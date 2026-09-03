import argparse
from datetime import datetime, timedelta
import netCDF4 as nc
import numpy as np

# Modified Julian Day epoch (days since 1858-11-17 00:00:00)
MJD_EPOCH = datetime(1858, 11, 17, 0, 0, 0)


def datetime_to_mjd(dt: datetime) -> float:
    return (dt - MJD_EPOCH).total_seconds() / 86400.0


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a multi-timestep FVCOM forcing file (given timestamp "
        "and hourly intervals), so bracket-based time interpolation works."
    )
    parser.add_argument(
        "--timestamp",
        type=str,
        default="2018-08-05T00:00:00.0000",
        help="First timestamp in ISO format (e.g., 2018-08-05T00:00:00.0000).",
    )
    parser.add_argument(
        "--timesteps",
        type=int,
        default=24,
        help="Number of hourly timesteps to generate (default: 24).",
    )
    parser.add_argument(
        "--source-file",
        type=str,
        default="gl_forcing_dummy.nc",
        help="Path to source NetCDF file for XLAT and XLONG.",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default="gl_fvcom_forcing_output.nc",
        help="Path for output NetCDF file.",
    )
    return parser.parse_args()


def parse_timestamp(ts_str: str) -> datetime:
    """Parse various ISO and standard date format strings into a datetime object."""
    formats = [
        "%Y-%m-%dT%H:%M:%S.%f",  # 2018-08-05T00:00:00.0000
        "%Y-%m-%dT%H:%M:%S",  # 2018-08-05T00:00:00
        "%Y-%m-%d %H:%M:%S",  # 2018-08-05 00:00:00
        "%Y-%m-%d",  # 2018-08-05
    ]
    for fmt in formats:
        try:
            return datetime.strptime(ts_str, fmt)
        except ValueError:
            continue
    raise ValueError(f"Unable to parse timestamp string: '{ts_str}'")


def format_times_str(dt: datetime) -> str:
    """Format a datetime as an ISO timestamp padded/truncated to exactly 26 chars."""
    return dt.strftime("%Y-%m-%dT%H:%M:%S.%f")[:26].ljust(26)


def main():
    args = parse_args()

    # Parse user-supplied timestamp and derive the sequence of hourly timesteps
    dt0 = parse_timestamp(args.timestamp)
    timesteps = [dt0 + timedelta(hours=i) for i in range(args.timesteps + 1)]

    mjd_values = [datetime_to_mjd(dt) for dt in timesteps]
    times_strs = [format_times_str(dt) for dt in timesteps]

    # Read XLAT, XLONG, and coordinates (x, y) if present from source file
    with nc.Dataset(args.source_file, "r") as src:
        xlat_data = src.variables["XLAT"][:]
        xlong_data = src.variables["XLONG"][:]
        n_south_north, n_west_east = xlat_data.shape
        x_data = (
            src.variables["x"][:]
            if "x" in src.variables
            else np.zeros((n_south_north, n_west_east), dtype=np.float32)
        )
        y_data = (
            src.variables["y"][:]
            if "y" in src.variables
            else np.zeros((n_south_north, n_west_east), dtype=np.float32)
        )

    # Initialize destination NetCDF file
    with nc.Dataset(args.output_file, "w", format="NETCDF4_CLASSIC") as dst:
        # Global Attributes
        dst.title = "FVCOM wrf_grid forcing"
        dst.source = (
            "wrf2fvcom version 0.13 (2007-07-19) (Bulk method: COARE 2.6Z)"
        )
        dst.modeler = "Xue's Group @MTU"
        dst.history = (
            f"Generated via Python netCDF script for {args.timesteps}-timestep forcing "
            "(bracket interpolation requires >= 2 time records)."
        )

        # Dimensions
        dst.createDimension("time", None)  # Unlimited dimension
        dst.createDimension("south_north", n_south_north)
        dst.createDimension("west_east", n_west_east)
        dst.createDimension("DateStrLen", 26)

        # Coordinate Variables (2D/1D)
        v_xlat = dst.createVariable(
            "XLAT", "f4", ("south_north", "west_east")
        )
        v_xlat.long_name = "latitude"
        v_xlat.description = "latitude,south is negative"
        v_xlat.units = "degrees_north"
        v_xlat[:] = xlat_data

        v_xlong = dst.createVariable(
            "XLONG", "f4", ("south_north", "west_east")
        )
        v_xlong.long_name = "longitude"
        v_xlong.description = "longitude,west is negative"
        v_xlong.units = "degrees_east"
        v_xlong[:] = xlong_data

        v_x = dst.createVariable("x", "f4", ("south_north", "west_east"))
        v_x.long_name = "nodal x-coordinate"
        v_x.units = "meters"
        v_x[:] = x_data

        v_y = dst.createVariable("y", "f4", ("south_north", "west_east"))
        v_y.long_name = "nodal y-coordinate"
        v_y.units = "meters"
        v_y[:] = y_data

        # Time Variables -- written for all timesteps

        v_time = dst.createVariable("time", "f4", ("Time",))
        v_time.long_name = "time"
        v_time.units = "days since 0.0"
        v_time.format = "relative time for forcing"
        v_time.time_zone = "none"

        for i in range(args.timesteps + 1):
            v_time[i] = i / 24.0

        # Data Variables (3D: Time, south_north, west_east)
        fields = [
            (
                "Longwave",
                "long wave radiation",
                "W/m2",
                "longwave,upward is negative",
            ),
            ("SPQ", "specific humidity", "kg/kg", "sea surface specific humidity"),
            ("Shortwave", "short wave radiation", "W/m2", "upward is negative"),
            ("U10", "eastward wind speed", "m/s", "U at 10m"),
            ("V10", "northward wind speed", "m/s", "V at 10m"),
            ("air_pressure", "air pressure", "Pa", "sea surface air pressure"),
            (
                "air_temperature",
                "air temperature",
                "degree (C)",
                "sea surface air temperature",
            ),
            ("cloud_cover", "cloud cover", "", "cloud cover,0-1"),
            ("dew_point", "dew point", "degree (C)", "dew point temperature"),
            (
                "relative_humidity",
                "relative humidity",
                "kg/kg",
                "sea surface relative humidity",
            ),
        ]
        
        # Create a full 3D array of zeros across the temporal dimension
        zero_grid_3d = np.zeros((args.timesteps + 1, n_south_north, n_west_east), dtype=np.float32)
        
        for name, long_name, units, desc in fields:
            var = dst.createVariable(
                name, "f4", ("Time", "south_north", "west_east")
            )
            var.long_name = long_name
            var.units = units
            var.description = desc
            var.coordinates = "XLONG XLAT"
            # Assign the entire 3D block at once
            var[:] = zero_grid_3d

    ts_summary = f"{times_strs[0].strip()} to {times_strs[-1].strip()}"
    print(f"Successfully created '{args.output_file}' with {args.timesteps + 1} timesteps: {ts_summary}")


if __name__ == "__main__":
    main()
