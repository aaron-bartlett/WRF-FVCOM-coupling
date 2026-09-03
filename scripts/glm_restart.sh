#!/bin/bash --login
#SBATCH -A glm200001
#SBATCH --job-name=glm_restart
#SBATCH --time=01:00:00
#SBATCH --output=/compass/glm200001/cmu/coupled-run/logs/%x_%j.out
#SBATCH --error=/compass/glm200001/cmu/coupled-run/logs/%x_%j.err

# glm_restart.sh -- top-level driver for the coupled WRF/FVCOM Great Lakes run.
#
# Meant to be invoked repeatedly (e.g. from cron).  Each invocation:
#   1. exits immediately if a coupled_run_glm / real_glm / metgrid_glm job is
#      already in the queue
#   2. derives the latest WRF and FVCOM restart timestamps; if they disagree it
#      logs that and moves forward with the earliest
#   3. runs check_wrf_inputs.sh:
#        exit 0  -> inputs still cover the window: re-arm WRF + FVCOM namelists
#                   for a restart and submit the coupled run
#        exit 10 -> check_wrf_inputs.sh launched the WPS pipeline; nothing else
#                   to do this cycle
#   4. always runs wget_cdsapi_requests.sh
#
# All output is appended to log.glm_restart.
set -uo pipefail

SCRIPT_DIR=/compass/glm200001/cmu/coupled-run/scripts
ROOT=/compass/glm200001/cmu/coupled-run
export GLM_COUPLED_ROOT="$ROOT"
LOG="$ROOT/log.glm_restart"
export GLM_LOG="$LOG"

WRF_RUN="$ROOT/nu-wrf-v11_cpl_oasis4/WRF/run"
FVCOM_RUN="$ROOT/FVCOM41_oasis_wrf_fvcom_iceDynamic_new/run"
WRF_NML="$WRF_RUN/namelist.input"
FVCOM_NML="$FVCOM_RUN/gl_run.nml"
COUPLED_SUBMIT="$ROOT/submit_coupledrun.sh"

SQUEUE_BIN="$(command -v squeue || true)"; SQUEUE_BIN="${SQUEUE_BIN:-squeue}"
SBATCH_BIN="$(command -v sbatch || true)"; SBATCH_BIN="${SBATCH_BIN:-sbatch}"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG" >&2; }
finish() { 
  log "Submitting wget_cdsapi_requests.sh via sbatch..." 
  sbatch "$SCRIPT_DIR/wget_cdsapi_requests.sh" >>"$LOG" 2>&1 || log "sbatch submission failed" 
  log "================ glm_restart.sh end ================" 
  exit "${1:-0}" 
}

log "================ glm_restart.sh start ================"

# ---- 1. don't stack jobs ---------------------------------------------
running="$("$SQUEUE_BIN" -u "${USER:-$(whoami)}" -h -o '%j' 2>/dev/null \
    | grep -Ewx 'coupled_run_glm|real_glm|metgrid_glm' || true)"
if [[ -n "$running" ]]; then
    log "job(s) already queued: $(echo "$running" | sort -u | tr '\n' ' ')- exiting, no action"
    log "================ glm_restart.sh end (queue busy) ================"
    exit 0
fi

command -v ncdump >/dev/null 2>&1 || source "$ROOT/load_modules.sh" >>"$LOG" 2>&1 || true

# ---- 2. latest restart timestamps ----------------------------------
latest_wrfrst="$(ls -t "$WRF_RUN"/wrfrst_d01_* 2>/dev/null | head -1 || true)"
if [[ -z "$latest_wrfrst" ]]; then
    log "no wrfrst_d01_* in $WRF_RUN -- cannot determine WRF restart time"
    finish 1
fi
wrf_ts="$(basename "$latest_wrfrst" | sed -E 's/^wrfrst_d01_//')"
wrf_epoch="$(date -d "${wrf_ts//[T_]/ }" +%s 2>/dev/null || true)"
[[ -n "$wrf_epoch" ]] || { log "cannot parse WRF restart ts from '$latest_wrfrst'"; finish 1; }
log "latest WRF   restart: $wrf_ts   ($(basename "$latest_wrfrst"))"

latest_fvrst="$(ls -t "$FVCOM_RUN"/output/gl_restart_*.nc 2>/dev/null | head -1 || true)"
if [[ -z "$latest_fvrst" ]]; then
    log "no output/gl_restart_*.nc in $FVCOM_RUN -- cannot determine FVCOM restart time"
    finish 1
fi
fv_raw="$(ncdump -v Times "$latest_fvrst" 2>/dev/null | tail -n 2 | head -n 1 \
    | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}' | head -1)"
[[ -n "$fv_raw" ]] || { log "cannot read Times from $latest_fvrst"; finish 1; }
fv_ts="${fv_raw/T/_}"
fv_epoch="$(date -d "${fv_ts//[T_]/ }" +%s 2>/dev/null || true)"
[[ -n "$fv_epoch" ]] || { log "cannot parse FVCOM restart ts '$fv_raw'"; finish 1; }
log "latest FVCOM restart: $fv_ts   ($(basename "$latest_fvrst"))"

fvrst_input="$FVCOM_RUN/input/$(basename "$latest_fvrst")"
if cp -f "$latest_fvrst" "$fvrst_input"; then
    log "staged FVCOM restart: $latest_fvrst -> $fvrst_input"
else
    log "failed to copy $latest_fvrst -> $fvrst_input"
    finish 1
fi


if (( wrf_epoch == fv_epoch )); then
    restart_epoch=$wrf_epoch
else
    log "WRF and FVCOM restart timestamps DIFFER (WRF $wrf_ts vs FVCOM $fv_ts); moving forward with earliest"
    restart_epoch=$(( wrf_epoch < fv_epoch ? wrf_epoch : fv_epoch ))
fi
RESTART_TS="$(date -d "@${restart_epoch}" +%Y-%m-%d_%H:%M:%S)"
log "restart timestamp: $RESTART_TS"

# ---- 3. validate WRF inputs / maybe launch WPS -------------------
"$SCRIPT_DIR/check_wrf_inputs.sh" "$RESTART_TS"
cwi_rc=$?

case "$cwi_rc" in
  0)
    log "WRF inputs valid -- re-arming namelists for restart and submitting coupled run"

    # coupling window length currently configured in namelist.input
    _num() { grep -iE "^[[:space:]]*$1[[:space:]]*=" "$WRF_NML" | grep -oE '[0-9]+' | head -1; }
    rd="$(_num run_days)";    rd="${rd:-0}"
    rh="$(_num run_hours)";   rh="${rh:-0}"
    rmn="$(_num run_minutes)"; rmn="${rmn:-0}"
    rsc="$(_num run_seconds)"; rsc="${rsc:-0}"
    dur=$(( rd*86400 + rh*3600 + rmn*60 + rsc ))
    (( dur > 0 )) || dur=43200          # fall back to namcouple $RUNTIME (12 h)
    end_epoch=$(( restart_epoch + dur ))

    wrf_set() {
        grep -qiE "^[[:space:]]*$1[[:space:]]*=" "$WRF_NML" || { log "WRF namelist: key $1 missing"; return 1; }
        sed -i -E "s|^([[:space:]]*$1[[:space:]]*=[[:space:]]*).*|\1$2|I" "$WRF_NML"
    }
    wrf_set start_year  "$(date -d "@$restart_epoch" +%Y),"
    wrf_set start_month "$(date -d "@$restart_epoch" +%m),"
    wrf_set start_day   "$(date -d "@$restart_epoch" +%d),"
    wrf_set start_hour  "$(date -d "@$restart_epoch" +%H),"
    wrf_set end_year    "$(date -d "@$end_epoch" +%Y),"
    wrf_set end_month   "$(date -d "@$end_epoch" +%m),"
    wrf_set end_day     "$(date -d "@$end_epoch" +%d),"
    wrf_set end_hour    "$(date -d "@$end_epoch" +%H),"
    wrf_set run_days    "0,"
    wrf_set run_hours   "$(( dur / 3600 )),"
    wrf_set run_minutes "$(( (dur % 3600) / 60 )),"
    wrf_set run_seconds "$(( dur % 60 )),"
    wrf_set restart     ".true.,"

    # FVCOM restart re-arm (keep both coupled components consistent)
    fv_base="$(basename "$latest_fvrst")"
    fv_start="$(date -d "@$restart_epoch" '+%Y-%m-%d %H:%M:%S')"
    fv_end="$(date -d "@$end_epoch" '+%Y-%m-%d %H:%M:%S')"
    fv_set() { sed -i -E "s|^([[:space:]]*$1[[:space:]]*=[[:space:]]*).*|\1$2|I" "$FVCOM_NML"; }
    fv_set START_DATE     "'${fv_start}',"
    fv_set END_DATE       "'${fv_end}'"
    fv_set STARTUP_TYPE   "'hotstart',"
    fv_set STARTUP_FILE   "'./${fv_base}',"
    fv_set RST_FIRST_OUT  "'${fv_start}',"
    fv_set NC_FIRST_OUT   "'${fv_start}',"
    fv_set NCAV_FIRST_OUT "'${fv_start}',"
    log "re-armed: WRF $(date -d "@$restart_epoch" +%Y-%m-%d_%H:%M:%S) -> $(date -d "@$end_epoch" +%Y-%m-%d_%H:%M:%S) restart=.true.; FVCOM hotstart './output/${fv_base}'"

    if [[ -f "$COUPLED_SUBMIT" ]]; then
        if jid="$("$SBATCH_BIN" --parsable "$COUPLED_SUBMIT" 2>>"$LOG")"; then
            log "submitted coupled run job ${jid%%;*} (job-name coupled_run_glm)"
        else
            log "ERROR: sbatch $COUPLED_SUBMIT failed"
        fi
    else
        log "ERROR: $COUPLED_SUBMIT not found -- coupled run not submitted"
    fi
    ;;
  10)
    log "check_wrf_inputs.sh launched the WPS pipeline; coupled run deferred to a later cycle"
    ;;
  *)
    log "check_wrf_inputs.sh exited $cwi_rc (unexpected); coupled run not submitted"
    ;;
esac

# ---- 4. always: pull any ready CDS downloads --------------------
finish 0
