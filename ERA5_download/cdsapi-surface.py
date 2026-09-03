#!/usr/bin/env python3
"""Submit ONE asynchronous CDS request for a full year of ERA5 single-level data.

Called once per year by submit_ERA5_download.sh.  The request is submitted with
wait_until_complete=False so the CDS request_id is returned immediately; the
actual file is pulled later by wget_cdsapi_requests.sh.

Machine-readable output (parsed by the caller):
    REQUEST_ID=<id>
    TARGET=<filename>
"""
from pathlib import Path
import argparse
import sys

import cdsapi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True)
    args = ap.parse_args()

    year = args.year
    target = f"era5_surf_{year}.grib"

    if Path(target).is_file():
        print(f"TARGET={target}")
        print("REQUEST_ID=EXISTS")
        return

    dataset = "reanalysis-era5-single-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "2m_dewpoint_temperature",
            "2m_temperature",
            "mean_sea_level_pressure",
            "sea_surface_temperature",
            "surface_pressure",
            "total_precipitation",
            "skin_temperature",
            "surface_latent_heat_flux",
            "top_net_solar_radiation_clear_sky",
            "snow_depth",
            "soil_temperature_level_1",
            "soil_temperature_level_2",
            "soil_temperature_level_3",
            "soil_temperature_level_4",
            "soil_type",
            "volumetric_soil_water_layer_1",
            "volumetric_soil_water_layer_2",
            "volumetric_soil_water_layer_3",
            "volumetric_soil_water_layer_4",
            "leaf_area_index_high_vegetation",
            "geopotential",
            "land_sea_mask",
            "sea_ice_cover",
        ],
        "year": [f"{year}"],
        "month": ["01", "02", "03", "04", "05", "06",
                  "07", "08", "09", "10", "11", "12"],
        "day": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"],
        "time": ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"],
        "data_format": "grib",
        "download_format": "unarchived",
        "area": [55, -103, 35, -67],  # North, West, South, East
    }

    c = cdsapi.Client(wait_until_complete=False, delete=False)
    result = c.retrieve(dataset, request)
    req_id = result.reply.get("request_id")

    print(f"TARGET={target}")
    print(f"REQUEST_ID={req_id}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001 - surface any failure to the caller
        print("REQUEST_ID=ERROR")
        print(f"# cdsapi-surface.py failed: {exc}", file=sys.stderr)
        sys.exit(1)
