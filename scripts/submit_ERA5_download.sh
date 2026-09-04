#!/bin/bash --login
# submit_ERA5_download.sh --input-year YYYY [--input-month MM]
#
# Fire asynchronous CDS requests for ERA5 forcing:
#   * pressure-level data : one request per month -- every month 1..12, or just
#     --input-month MM when given (up to 12 requests, or exactly 1)
#   * single-level  data  : one request for the whole YYYY year (1 request)
#
# --input-month is used by glm_restart.sh's proactive look-ahead request (fires
# every cycle for the single month/year 5 years past the current restart time,
# rather than pulling the whole year at once); submit_WPS.sh's own reactive
# calls omit it and get the full 12-month year, as before.
#
# Each newly submitted request is appended to ERA5_download/cdsapi_requests.csv.
# The requests use wait_until_complete=False, so nothing is downloaded here --
# wget_cdsapi_requests.sh pulls the finished files on a later cycle.
#
# A request is skipped when either the target grib already exists on disk or an
# undownloaded row for the same year/month is already in the csv.
set -uo pipefail

SCRIPT_DIR=/compass/glm200001/cmu/coupled-run/scripts
ROOT=/compass/glm200001/cmu/coupled-run
ERA5_DIR="$ROOT/ERA5_download"
LOG="${GLM_LOG:-$ROOT/log.glm_restart}"
CSV="$ERA5_DIR/cdsapi_requests.csv"
PLEVS_DIR="$ERA5_DIR/plevs-ERA5"
SURF_DIR="$ERA5_DIR/surface-ERA5"

module load py-pip
source /home/bart753/.venv/bin/activate
PYTHON_BIN=python3

log() { printf '[%s] submit_ERA5_download: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG" >&2; }

YEAR=""
MONTH=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-year)    YEAR="${2:-}"; shift 2;;
        --input-year=*)  YEAR="${1#*=}"; shift;;
        --input-month)   MONTH="${2:-}"; shift 2;;
        --input-month=*) MONTH="${1#*=}"; shift;;
        -h|--help)       echo "usage: $0 --input-year YYYY [--input-month MM]"; exit 0;;
        *)               log "unknown argument: $1"; exit 2;;
    esac
done
[[ "$YEAR" =~ ^[0-9]{4}$ ]] || { log "usage: $0 --input-year YYYY [--input-month MM]"; exit 2; }
if [[ -n "$MONTH" ]]; then
    [[ "$MONTH" =~ ^(0?[1-9]|1[0-2])$ ]] \
        || { log "usage: $0 --input-year YYYY [--input-month MM] (MM must be 1-12)"; exit 2; }
fi

mkdir -p "$PLEVS_DIR/plevs-ERA5-$YEAR" "$SURF_DIR"
[[ -f "$CSV" ]] || echo "request_timestamp,request_year,request_month,request_id,downloaded" > "$CSV"

# true (exit 0) when an undownloaded request row already exists for year/month.
pending_row() {  # year month(may be "")
    awk -F, -v y="$1" -v m="$2" '
        NR > 1 && $2 == y && $3 == m && $5 == 0 &&
        $4 != "" && $4 != "ERROR" && $4 != "EXISTS" { found = 1 }
        END { exit !found }' "$CSV"
}

append_row() {  # year month request_id
    printf '%s,%s,%s,%s,0\n' "$(date '+%Y-%m-%dT%H:%M:%S')" "$1" "$2" "$3" >> "$CSV"
}

# Submit one request; echoes the request_id on success, nothing otherwise.
submit_one() {  # label pyscript arg...
    local label="$1"; shift
    local out rid rc
    out="$(cd "$ERA5_DIR" && "$PYTHON_BIN" "$@" 2>>"$LOG")"; rc=$?
    rid="$(printf '%s\n' "$out" | sed -n 's/^REQUEST_ID=//p' | tail -1)"
    if [[ $rc -ne 0 || -z "$rid" || "$rid" == "ERROR" ]]; then
        log "$label: request submission FAILED"
        return 1
    fi
    if [[ "$rid" == "EXISTS" ]]; then
        log "$label: target already present, skipping"
        return 2
    fi
    printf '%s' "$rid"
}

# ---- pressure levels: one request per month (or just --input-month) --------
if [[ -n "$MONTH" ]]; then
    month_list=("$((10#$MONTH))")
else
    month_list=($(seq 1 12))
fi
for m in "${month_list[@]}"; do
    mm="$(printf '%02d' "$m")"
    target="$PLEVS_DIR/plevs-ERA5-$YEAR/era5_plevs_${YEAR}-${mm}.grib"
    if [[ -f "$target" ]]; then
        log "plevs $YEAR-$mm: file exists, skipping"
        continue
    fi
    if pending_row "$YEAR" "$mm"; then
        log "plevs $YEAR-$mm: undownloaded request already in csv, skipping"
        continue
    fi
    if rid="$(submit_one "plevs $YEAR-$mm" cdsapi-levels.py --year "$YEAR" --month "$m")"; then
        append_row "$YEAR" "$mm" "$rid"
        log "plevs $YEAR-$mm: submitted request_id=$rid"
    fi
done

# ---- single levels: one request for the whole year ------------------------
target="$SURF_DIR/era5_surf_${YEAR}.grib"
if [[ -f "$target" ]]; then
    log "surface $YEAR: file exists, skipping"
elif pending_row "$YEAR" ""; then
    log "surface $YEAR: undownloaded request already in csv, skipping"
elif rid="$(submit_one "surface $YEAR" cdsapi-surface.py --year "$YEAR")"; then
    append_row "$YEAR" "" "$rid"
    log "surface $YEAR: submitted request_id=$rid"
fi

log "done for year $YEAR${MONTH:+ (month $MONTH only for pressure levels)}"
