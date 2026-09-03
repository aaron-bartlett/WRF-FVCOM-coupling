#!/usr/bin/env python3
"""Submit ONE asynchronous CDS request for ERA5 pressure-level data.

Called once per (year, month) by submit_ERA5_download.sh.  The request is
submitted with wait_until_complete=False so the CDS request_id is returned
immediately; the actual file is pulled later by wget_cdsapi_requests.sh.

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
    ap.add_argument("--month", type=int, required=True)
    args = ap.parse_args()

    year = args.year
    month = args.month
    target = f"era5_plevs_{year}-{month:02d}.grib"

    if Path(target).is_file():
        print(f"TARGET={target}")
        print("REQUEST_ID=EXISTS")
        return

    dataset = "reanalysis-era5-pressure-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": [
            "geopotential",
            "relative_humidity",
            "specific_humidity",
            "temperature",
            "u_component_of_wind",
            "v_component_of_wind",
        ],
        "year": [f"{year}"],
        "month": [f"{month:02d}"],
        "day": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"],
        "time": ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"],
        "pressure_level": [
            "1", "2", "3", "5", "7", "10", "20", "30", "50", "70", "100", "125",
            "150", "175", "200", "225", "250", "300", "350", "400", "450", "500",
            "550", "600", "650", "700", "750", "775", "800", "825", "850", "875",
            "900", "925", "950", "975", "1000",
        ],
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
        print(f"REQUEST_ID=ERROR", file=sys.stdout)
        print(f"# cdsapi-levels.py failed: {exc}", file=sys.stderr)
        sys.exit(1)
