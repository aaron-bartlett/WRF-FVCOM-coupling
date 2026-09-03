#!/bin/bash --login
#SBATCH -A glm200001
#SBATCH --job-name=glm_cdsapi_dl
#SBATCH --time=01:00:00
#SBATCH --output=/compass/glm200001/cmu/coupled-run/logs/%x_%j.out
#SBATCH --error=/compass/glm200001/cmu/coupled-run/logs/%x_%j.err


# wget_cdsapi_requests.sh
#
# Walk ERA5_download/cdsapi_requests.csv.  For every row still flagged
# not-downloaded (downloaded == 0), ask the CDS API whether the job has
# finished and, if so, wget the result into the correct ERA5_download
# sub-directory and flip the flag to 1.
#
# Data that is missing / not yet ready is skipped silently -- no error, no
# warning.  Successful downloads are recorded in log.cdsapi_downloads.
set -uo pipefail

SCRIPT_DIR=/compass/glm200001/cmu/coupled-run/scripts
ROOT=/compass/glm200001/cmu/coupled-run
ERA5_DIR="$ROOT/ERA5_download"
CSV="$ERA5_DIR/cdsapi_requests.csv"
DL_LOG="$ROOT/log.cdsapi_downloads"
RCFILE="${CDSAPI_RC:-$HOME/.cdsapirc}"

[[ -f "$CSV" ]] || exit 0
[[ -f "$RCFILE" ]] || { echo "wget_cdsapi_requests: $RCFILE not found" >&2; exit 0; }

CDS_URL="$(sed -n 's/^[[:space:]]*url[[:space:]]*:[[:space:]]*//p' "$RCFILE" | head -1)"
CDS_KEY="$(sed -n 's/^[[:space:]]*key[[:space:]]*:[[:space:]]*//p' "$RCFILE" | head -1)"
CDS_URL="${CDS_URL%/}"
CDS_URL="${CDS_URL%%[[:space:]]}"
CDS_KEY="${CDS_KEY%%[[:space:]]}"
[[ -n "$CDS_URL" && -n "$CDS_KEY" ]] || { echo "wget_cdsapi_requests: could not parse url/key from $RCFILE" >&2; exit 0; }

dl_log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >> "$DL_LOG"; }

# Echo the download href for a finished job, or nothing.
job_href() {  # request_id
    local rid="$1" json
    json="$(wget -q -O- --header="PRIVATE-TOKEN: ${CDS_KEY}" \
            "${CDS_URL}/retrieve/v1/jobs/${rid}/results" 2>/dev/null)" || return 0
    [[ -n "$json" ]] || return 0
    printf '%s' "$json" | python3 -c '
import sys, json
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
if not isinstance(d, dict):
    sys.exit(0)
v = d.get("asset", {}).get("value", {})
print(v.get("href", "") or "")
' 2>/dev/null
}

TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT
changed=0
first=1

while IFS=, read -r ts year month rid flag || [[ -n "${ts:-}" ]]; do
    # header row / blank line -> copy through untouched
    if [[ $first -eq 1 ]]; then
        first=0
        printf '%s,%s,%s,%s,%s\n' "$ts" "$year" "$month" "$rid" "$flag" >> "$TMP"
        continue
    fi
    if [[ -z "${ts:-}${rid:-}" ]]; then
        continue
    fi

    if [[ "${flag:-0}" != "0" || -z "${rid:-}" || "$rid" == "ERROR" || "$rid" == "EXISTS" ]]; then
        printf '%s,%s,%s,%s,%s\n' "$ts" "$year" "$month" "$rid" "$flag" >> "$TMP"
        continue
    fi

    # Destination is fully determined by the rest of the row.
    if [[ -n "${month// /}" ]]; then
        dest="$ERA5_DIR/plevs-ERA5/plevs-ERA5-${year}/era5_plevs_${year}-${month}.grib"
    else
        dest="$ERA5_DIR/surface-ERA5/era5_surf_${year}.grib"
    fi

    if [[ -f "$dest" ]]; then
        printf '%s,%s,%s,%s,1\n' "$ts" "$year" "$month" "$rid" >> "$TMP"
        changed=1
        continue
    fi

    href="$(job_href "$rid")"
    if [[ -z "$href" ]]; then
        # not ready / not available -- leave the row alone, stay quiet
        printf '%s,%s,%s,%s,0\n' "$ts" "$year" "$month" "$rid" >> "$TMP"
        continue
    fi

    mkdir -p "$(dirname "$dest")"
    if wget -q --header="PRIVATE-TOKEN: ${CDS_KEY}" -O "${dest}.part" "$href"; then
        mv "${dest}.part" "$dest"
        printf '%s,%s,%s,%s,1\n' "$ts" "$year" "$month" "$rid" >> "$TMP"
        changed=1
        dl_log "downloaded $(basename "$dest")  <- request $rid"
    else
        rm -f "${dest}.part"
        printf '%s,%s,%s,%s,0\n' "$ts" "$year" "$month" "$rid" >> "$TMP"
    fi
done < "$CSV"

if [[ "$changed" -eq 1 ]]; then
    cat "$TMP" > "$CSV"
fi
