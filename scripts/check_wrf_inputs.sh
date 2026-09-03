#!/bin/bash --login
# check_wrf_inputs.sh <restart_ts>
#   restart_ts : YYYY-MM-DD_HH:MM:SS  (timestamp the coupled run will restart from)
#
# Reads wrfbdy_d01's Times.  If the restart time is at or beyond the
# second-to-last boundary time, the current wrfbdy no longer covers the upcoming
# window -> run submit_WPS.sh and exit 10.
# Otherwise the existing WRF inputs still cover the window -> exit 0.
#
# exit codes: 0 = inputs valid, no WPS needed
#            10 = WPS pipeline launched
#             2 = usage / parse error
set -uo pipefail

SCRIPT_DIR=/compass/glm200001/cmu/coupled-run/scripts
ROOT=/compass/glm200001/cmu/coupled-run
LOG="${GLM_LOG:-$ROOT/log.glm_restart}"
WRF_RUN="$ROOT/nu-wrf-v11_cpl_oasis4/WRF/run"
WRFBDY="$WRF_RUN/wrfbdy_d01"

log() { printf '[%s] check_wrf_inputs: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG" >&2; }

RESTART_TS="${1:-}"
[[ -n "$RESTART_TS" ]] || { log "usage: $0 <YYYY-MM-DD_HH:MM:SS>"; exit 2; }

command -v ncdump >/dev/null 2>&1 || source "$ROOT/load_modules.sh" >>"$LOG" 2>&1 || true

restart_epoch="$(date -d "${RESTART_TS//[T_]/ }" +%s 2>/dev/null || true)"
[[ -n "$restart_epoch" ]] || { log "cannot parse restart ts '$RESTART_TS'"; exit 2; }

launch_wps() {
    log "launching WPS pipeline for restart $RESTART_TS"
    "$SCRIPT_DIR/submit_WPS.sh" "$RESTART_TS" || log "submit_WPS.sh returned non-zero"
    exit 10
}

if [[ ! -f "$WRFBDY" ]]; then
    log "wrfbdy_d01 not found -> WPS required"
    launch_wps
fi

# every boundary timestamp, in file order
BDY_TIMES=()
while IFS= read -r _t; do
    [[ -n "$_t" ]] && BDY_TIMES+=("$_t")
done < <(
    ncdump -v Times "$WRFBDY" 2>/dev/null \
        | sed -n '/Times[[:space:]]*=/,/;/p' \
        | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}[T_][0-9]{2}:[0-9]{2}:[0-9]{2}'
)
n=${#BDY_TIMES[@]}
if (( n < 2 )); then
    log "read $n (<2) Times entries from wrfbdy_d01 -> WPS required"
    launch_wps
fi

penult="${BDY_TIMES[n-2]}"
last="${BDY_TIMES[n-1]}"
penult_epoch="$(date -d "${penult//[T_]/ }" +%s 2>/dev/null || true)"
[[ -n "$penult_epoch" ]] || { log "cannot parse wrfbdy time '$penult'"; exit 2; }

log "restart=$RESTART_TS  wrfbdy covers .. $penult (second-to-last) / $last (last)"

if (( restart_epoch >= penult_epoch )); then
    log "restart is at/after the second-to-last boundary time -> WPS required"
    launch_wps
fi

log "existing WRF inputs still valid; no WPS needed"
exit 0
