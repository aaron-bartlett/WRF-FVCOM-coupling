#!/bin/bash --login
# rename_fvcom_output.sh
#
# FVCOM writes its NetCDF output and restart files with an opaque sequential
# suffix (gl_0001.nc, gl_0002.nc, ... and gl_restart_0001.nc, ...). This scans
# FVCOM_RUN/output and renames every such file to a timestamped name derived
# from the FIRST entry of that file's own Times variable:
#
#   gl_0001.nc         -> gl_2024-05-01_00:00:00.nc
#   gl_restart_0001.nc -> gl_restart_2024-05-01_00:00:00.nc
#
# so files identify and sort by simulation time instead of an opaque sequence
# number. Run once per glm_restart.sh cycle, before the latest restart file is
# located.
#
# Usage: ./rename_fvcom_output.sh
#
# Idempotent: files that are already timestamped, or don't match the
# sequential gl_NNNN.nc / gl_restart_NNNN.nc naming, are left alone. A rename
# that would overwrite an existing file is skipped with a warning instead.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${GLM_COUPLED_ROOT:-$(dirname "$SCRIPT_DIR")}"
LOG="${GLM_LOG:-$ROOT/log.glm_restart}"
FVCOM_RUN="$ROOT/FVCOM41_oasis_wrf_fvcom_iceDynamic_new/run"
OUT_DIR="$FVCOM_RUN/output"

log() { printf '[%s] rename_fvcom_output: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG" >&2; }

command -v ncdump >/dev/null 2>&1 || source "$ROOT/load_modules.sh" >>"$LOG" 2>&1 || true

# rename_one <file> : rename <file> to a timestamped name based on the first
# entry of its own Times variable. Preserves the gl_ / gl_restart_ prefix.
rename_one() {
    local f="$1" base dir first_raw first_ts new dest
    base="$(basename "$f")"
    dir="$(dirname "$f")"

    first_raw="$(ncdump -v Times "$f" 2>/dev/null \
        | sed -n '/Times[[:space:]]*=/,/;/p' \
        | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}[T_ ][0-9]{2}:[0-9]{2}:[0-9]{2}' \
        | head -1)"
    if [[ -z "$first_raw" ]]; then
        log "WARNING: could not read Times from $f -- leaving as-is"
        return 1
    fi
    first_ts="$(date -d "${first_raw//[T_]/ }" +%Y-%m-%d_%H:00:00 2>/dev/null || true)"
    if [[ -z "$first_ts" ]]; then
        log "WARNING: could not parse Times '$first_raw' from $f -- leaving as-is"
        return 1
    fi

    case "$base" in
        gl_restart_*) new="gl_restart_${first_ts}.nc" ;;
        *)            new="gl_${first_ts}.nc" ;;
    esac
    [[ "$base" == "$new" ]] && return 0   # already named correctly

    dest="$dir/$new"
    if [[ -e "$dest" ]]; then
        log "WARNING: $dest already exists -- not overwriting, leaving $base as-is"
        return 1
    fi
    if mv -n -- "$f" "$dest"; then
        log "renamed $base -> $new"
    else
        log "ERROR: failed to rename $f -> $dest"
        return 1
    fi
}

if [[ ! -d "$OUT_DIR" ]]; then
    log "no output dir $OUT_DIR -- nothing to rename"
    exit 0
fi

shopt -s nullglob
for f in "$OUT_DIR"/gl_*.nc; do
    base="$(basename "$f")"
    # only touch the opaque sequential names; already-timestamped or
    # otherwise-named files are left untouched.
    if [[ "$base" =~ ^gl_restart_[0-9]+\.nc$ ]] || [[ "$base" =~ ^gl_[0-9]+\.nc$ ]]; then
        rename_one "$f"
    fi
done
shopt -u nullglob

exit 0
