import argparse
from datetime import datetime
import netCDF4 as nc
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a multi-timestep FVCOM forcing file with time, Itime, and Itime2."
    )
    parser.add_argument(
        "--start-time",
        type=str,
        default="1996-01-01T00:00:00",
        help="Reference timestamp in ISO format (e.g., 1997-01-01T00:00:00).",
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
        default="gl_fvcom_forcing_dummy.nc",
        help="Path for output NetCDF file.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Generating 2 entries of periodic forcing, referenced to: {args.start_time}")

    # Read XLAT, XLONG, and coordinates (x, y) if present from source file
    try:
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
    except FileNotFoundError:
        print(f"Warning: Source file '{args.source_file}' not found. Using minimal dummy grid (10x10).")
        n_south_north, n_west_east = 10, 10
        xlat_data = np.zeros((n_south_north, n_west_east), dtype=np.float32)
        xlong_data = np.zeros((n_south_north, n_west_east), dtype=np.float32)
        x_data = np.zeros((n_south_north, n_west_east), dtype=np.float32)
        y_data = np.zeros((n_south_north, n_west_east), dtype=np.float32)

    # Initialize destination NetCDF file
    with nc.Dataset(args.output_file, "w", format="NETCDF4_CLASSIC") as dst:
        # Global Attributes
        dst.title = "FVCOM periodic forcing"
        dst.source = "wrf2fvcom version"
        dst.history = f"Generated for dummy periodic forcing referenced to {args.start_time}."

        # Dimensions
        dst.createDimension("time", None)  # Lowercase, unlimited
        dst.createDimension("south_north", n_south_north)
        dst.createDimension("west_east", n_west_east)

        # Coordinate Variables (2D)
        v_xlat = dst.createVariable("XLAT", "f4", ("south_north", "west_east"))
        v_xlat.long_name = "latitude"
        v_xlat.units = "degrees_north"
        v_xlat[:] = xlat_data

        v_xlong = dst.createVariable("XLONG", "f4", ("south_north", "west_east"))
        v_xlong.long_name = "longitude"
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

        # Time Variables (lowercase time dimension, explicit Itime/Itime2)
        v_time = dst.createVariable("time", "f4", ("time",))
        v_time.long_name = "time"
        v_time.units = "days since 0.0"
        v_time.time_zone = "none"
        v_time.format = "relative time for periodic forcing"

        v_itime = dst.createVariable("Itime", "i4", ("time",))
        v_itime.units = "days since 0.0"
        v_itime.time_zone = "none"
        v_itime.format = "relative integer days"

        v_itime2 = dst.createVariable("Itime2", "i4", ("time",))
        v_itime2.time_zone = "none"
        v_itime2.units = "msec since 00:00:00"

        v_time[0] = 0.0
        v_time[1] = 8766
        v_itime[0] = 0
        v_itime[1] = 8766
        v_itime2[0] = 0
        v_itime2[1] = 0

        # Data Variables (3D: time, south_north, west_east)
        fields = [
            ("Longwave", "long wave radiation", "W/m2", "longwave,upward is negative"),
            ("SPQ", "specific humidity", "kg/kg", "sea surface specific humidity"),
            ("Shortwave", "short wave radiation", "W/m2", "upward is negative"),
            ("U10", "eastward wind speed", "m/s", "U at 10m"),
            ("V10", "northward wind speed", "m/s", "V at 10m"),
            ("air_pressure", "air pressure", "Pa", "sea surface air pressure"),
            ("air_temperature", "air temperature", "degree (C)", "sea surface air temperature"),
            ("cloud_cover", "cloud cover", "", "cloud cover,0-1"),
            ("dew_point", "dew point", "degree (C)", "dew point temperature"),
            ("relative_humidity", "relative humidity", "kg/kg", "sea surface relative humidity"),
        ]
        
        # Create a full 3D array of zeros across the temporal dimension
        zero_grid_3d = np.zeros((2, n_south_north, n_west_east), dtype=np.float32)
        
        for name, long_name, units, desc in fields:
            var = dst.createVariable(name, "f4", ("time", "south_north", "west_east"))
            var.long_name = long_name
            var.units = units
            var.description = desc
            var.coordinates = "XLONG XLAT"
            var[:] = zero_grid_3d

    print(f"Successfully created '{args.output_file}' with 2 timesteps starting at relative time 0.0.")


if __name__ == "__main__":
    main()
