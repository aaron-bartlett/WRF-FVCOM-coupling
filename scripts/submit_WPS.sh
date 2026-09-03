#!/bin/bash --login
# submit_WPS.sh <restart_ts>
#   restart_ts : YYYY-MM-DD_HH:MM:SS
#
# Rebuild the WRF initial/boundary conditions for a fresh WPS window that
# starts at <restart_ts>:
#   1. verify every ERA5 grib file for the window is on disk
#      (else fire submit_ERA5_download.sh for each year in the window and bail)
#   2. ungrib the pressure-level gribs   (namelist.wps prefix = 'ERA5-PRES')
#   3. ungrib the single-level  gribs    (namelist.wps prefix = 'ERA5-SURF')
#   4. sbatch submit_metgrid.sh                    -> met_em files for the window
#   5. sbatch submit_real.sh (--dependency=afterok) -> wrfinput_d01 / wrfbdy_d01
#      (submit_real.sh refreshes the met_em links in WRF/run before real.exe)
#   6. widen WRF namelist.input + FVCOM gl_run.nml END_DATE to the new end date
#
# Staged / cron-driven: this returns once the batch jobs are queued.  While
# metgrid_glm / real_glm are in the queue, glm_restart.sh bails on its squeue
# guard; once real_glm finishes and wrfbdy_d01 covers the window, a later
# glm_restart.sh run submits the coupled run.
#
# * each ungrib step verifies the previous step succeeded before continuing.
set -uo pipefail

# ---- user-tunable settings --------------------------------------------------
# Length of the WPS window to build, in months, measured from <restart_ts>.
# (Was hard-coded at 4 years; 48 months is the equivalent default.)
WPS_WINDOW_MONTHS=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${GLM_COUPLED_ROOT:-$(dirname "$SCRIPT_DIR")}"
LOG="${GLM_LOG:-$ROOT/log.glm_restart}"
WPS_DIR="$ROOT/nu-wrf-v11_cpl_oasis4/WPS"
WRF_RUN="$ROOT/nu-wrf-v11_cpl_oasis4/WRF/run"
FVCOM_RUN="$ROOT/FVCOM41_oasis_wrf_fvcom_iceDynamic_new/run"
ERA5_DIR="$ROOT/ERA5_download"
WPS_NML="$WPS_DIR/namelist.wps"
WRF_NML="$WRF_RUN/namelist.input"
FVCOM_NML="$FVCOM_RUN/gl_run.nml"
METGRID_SUBMIT="$WPS_DIR/submit_metgrid.sh"
REAL_SUBMIT="$WRF_RUN/submit_real.sh"

# Verbose program output (module loads, link_grib.csh, ungrib.exe console) is
# written here, NOT into $LOG, so log.glm_restart keeps only the high-level
# progress lines emitted by log().  ungrib.exe also writes its own
# $WPS_DIR/ungrib.log, which stays the authoritative ungrib record.
WPS_LOG="$WPS_DIR/wps.log"

SBATCH_BIN="$(command -v sbatch || true)"; SBATCH_BIN="${SBATCH_BIN:-sbatch}"

log() { printf '[%s] submit_WPS: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG" >&2; }
die() { log "ERROR: $*"; exit 1; }

RESTART_TS="${1:-}"
[[ -n "$RESTART_TS" ]] || die "usage: $0 <YYYY-MM-DD_HH:MM:SS>"

start_h="${RESTART_TS//[T_]/ }"
date -d "$start_h" +%s >/dev/null 2>&1 || die "cannot parse restart ts '$RESTART_TS'"

# End of the window == start + WPS_WINDOW_MONTHS calendar months.
# NOTE: do the month arithmetic on the *date only* and re-attach the hour.
# GNU date parses a "+N" that immediately follows a HH:MM:SS time as a numeric
# timezone offset, not as a relative amount, so `date -d "$start_h +4 months"`
# silently means "9am at UTC+04:00, then +1 month" -- a ~1-month, tz-shifted
# window rather than a 4-month one. Splitting the time off avoids that.
_start_ymd="$(date -d "$start_h" +%Y-%m-%d)"
_start_hh="$( date -d "$start_h" +%H)"
end_h="$(date -d "${_start_ymd} +${WPS_WINDOW_MONTHS} months ${_start_hh}:00:00")"
date -d "$end_h" +%s >/dev/null 2>&1 || die "cannot compute window end from '$start_h' +$WPS_WINDOW_MONTHS months"

START_WPS="$(date -d "$start_h" +%Y-%m-%d_%H:00:00)"
END_WPS="$(  date -d "$end_h"   +%Y-%m-%d_%H:00:00)"
START_YR="$(date -d "$start_h" +%Y)"
END_YR="$(  date -d "$end_h"   +%Y)"
START_MO=$((10#$(date -d "$start_h" +%m)))
END_MO=$((  10#$(date -d "$end_h"   +%m)))
log "WPS window $START_WPS .. $END_WPS ($WPS_WINDOW_MONTHS months)"

# iterate (year month) pairs across the inclusive window; body uses $wy $wmm
for_each_month() {
    local wy=$START_YR wm=$START_MO wmm
    while (( wy < END_YR || (wy == END_YR && wm <= END_MO) )); do
        wmm="$(printf '%02d' "$wm")"
        "$@" "$wy" "$wmm"
        wm=$((wm + 1)); (( wm > 12 )) && { wm=1; wy=$((wy + 1)); }
    done
}

# ---- 1. ERA5 grib availability -------------------------------------------
missing=0
_check_plev() {  # year mm
    local f="$ERA5_DIR/plevs-ERA5/plevs-ERA5-$1/era5_plevs_$1-$2.grib"
    [[ -f "$f" ]] || { log "missing $f"; missing=1; }
}
for_each_month _check_plev
for (( yr = START_YR; yr <= END_YR; yr++ )); do
    f="$ERA5_DIR/surface-ERA5/era5_surf_${yr}.grib"
    [[ -f "$f" ]] || { log "missing $f"; missing=1; }
done

if (( missing )); then
    for (( yr = START_YR; yr <= END_YR; yr++ )); do
        log "requesting ERA5 download for $yr"
        "$SCRIPT_DIR/submit_ERA5_download.sh" --input-year "$yr" \
            || log "submit_ERA5_download.sh --input-year $yr returned non-zero"
    done
    die "ERA5 grib inputs incomplete for $START_WPS..$END_WPS; downloads requested, WPS aborted this cycle"
fi
log "all ERA5 gribs present for the window"

# ---- WPS environment ----------------------------------------------------
cd "$WPS_DIR" || die "cannot cd $WPS_DIR"
: >"$WPS_LOG" || true          # start this run's program-output log fresh
# shellcheck disable=SC1091
source "$ROOT/load_modules.sh" >>"$WPS_LOG" 2>&1 || die "load_modules.sh failed (see $WPS_LOG)"

# set a scalar key in a Fortran namelist: nml_set <file> <key> <value>
nml_set() {
    grep -qiE "^[[:space:]]*$2[[:space:]]*=" "$1" || die "key '$2' not found in $1"
    sed -i -E "s|^([[:space:]]*$2[[:space:]]*=[[:space:]]*).*|\1$3|I" "$1"
}

# run ungrib for one dataset: ungrib_dataset <prefix> <grib file...>
ungrib_dataset() {
    local prefix="$1"; shift
    nml_set "$WPS_NML" prefix "'${prefix}',"
    rm -f "$WPS_DIR"/GRIBFILE.* "$WPS_DIR/${prefix}:"* 2>/dev/null || true
    ./link_grib.csh "$@" >>"$WPS_LOG" 2>&1 || die "link_grib.csh ($prefix) failed"
    ./ungrib.exe >>"$WPS_LOG" 2>&1          || die "ungrib.exe ($prefix) failed (see $WPS_DIR/ungrib.log)"
    grep -q "Successful completion of program ungrib" "$WPS_DIR/ungrib.log" 2>/dev/null \
        || die "ungrib ($prefix) did not report success"
    ls "$WPS_DIR/${prefix}:"* >/dev/null 2>&1 \
        || die "ungrib ($prefix) produced no intermediate files"
    log "ungrib $prefix ok"
}

nml_set "$WPS_NML" start_date "'${START_WPS}',"
nml_set "$WPS_NML" end_date   "'${END_WPS}',"

# ---- 2. ungrib pressure levels ---------------------------------------
PLEVS_GRIB=()
_collect_plev() { PLEVS_GRIB+=("$ERA5_DIR/plevs-ERA5/plevs-ERA5-$1/era5_plevs_$1-$2.grib"); }
for_each_month _collect_plev
(( ${#PLEVS_GRIB[@]} > 0 )) || die "no pressure-level gribs resolved for the window"
ungrib_dataset "ERA5-PRES" "${PLEVS_GRIB[@]}"

# ---- 3. ungrib single levels ---------------------------------------
SURF_GRIB=()
for (( yr = START_YR; yr <= END_YR; yr++ )); do
    SURF_GRIB+=("$ERA5_DIR/surface-ERA5/era5_surf_${yr}.grib")
done
ungrib_dataset "ERA5-SURF" "${SURF_GRIB[@]}"

# ---- 4. metgrid (batch job: metgrid_glm) --------------------------
[[ -f "$METGRID_SUBMIT" ]] || die "$METGRID_SUBMIT not found"
MID="$("$SBATCH_BIN" --parsable "$METGRID_SUBMIT" "$START_WPS" "$END_WPS")" \
    || die "sbatch submit_metgrid.sh failed"
MID="${MID%%;*}"
log "queued metgrid job $MID (job-name metgrid_glm)"

# ---- 5. real (batch job: real_glm, after metgrid) ---------------
[[ -f "$REAL_SUBMIT" ]] || die "$REAL_SUBMIT not found"
RID="$("$SBATCH_BIN" --parsable --dependency=afterok:"$MID" "$REAL_SUBMIT" "$START_WPS" "$END_WPS")" \
    || die "sbatch submit_real.sh failed"
RID="${RID%%;*}"
log "queued real job $RID (job-name real_glm), depends on afterok:$MID"

# ---- 6. widen namelist windows to the new end date ---------------
s_y="$(date -d "$start_h" +%Y)"; s_m="$(date -d "$start_h" +%m)"
s_d="$(date -d "$start_h" +%d)"; s_h="$(date -d "$start_h" +%H)"
e_y="$(date -d "$end_h" +%Y)"; e_m="$(date -d "$end_h" +%m)"
e_d="$(date -d "$end_h" +%d)"; e_h="$(date -d "$end_h" +%H)"
win_days=$(( ( $(date -d "$end_h" +%s) - $(date -d "$start_h" +%s) ) / 86400 ))

nml_set "$WRF_NML" start_year  "${s_y},"
nml_set "$WRF_NML" start_month "${s_m},"
nml_set "$WRF_NML" start_day   "${s_d},"
nml_set "$WRF_NML" start_hour  "${s_h},"
nml_set "$WRF_NML" end_year    "${e_y},"
nml_set "$WRF_NML" end_month   "${e_m},"
nml_set "$WRF_NML" end_day     "${e_d},"
nml_set "$WRF_NML" end_hour    "${e_h},"
nml_set "$WRF_NML" run_days    "${win_days},"
nml_set "$WRF_NML" run_hours   "0,"
nml_set "$WRF_NML" run_minutes "0,"
nml_set "$WRF_NML" run_seconds "0,"
nml_set "$WRF_NML" restart     ".false.,"

# FVCOM gl_run.nml END_DATE == new WRF namelist end time
fvcom_end="$(date -d "$end_h" '+%Y-%m-%d %H:00:00')"
sed -i -E "s|^([[:space:]]*END_DATE[[:space:]]*=[[:space:]]*).*|\1'${fvcom_end}'|I" "$FVCOM_NML"
log "namelists widened: WRF end ${END_WPS} (${win_days}d), FVCOM END_DATE '${fvcom_end}'"

log "WPS pipeline queued; coupled run will be submitted by a later glm_restart.sh cycle"
exit 0
