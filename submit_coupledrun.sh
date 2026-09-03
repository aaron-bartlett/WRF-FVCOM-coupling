#!/bin/bash --login

#SBATCH -A glm200001
#SBATCH --job-name=coupled_run_glm
#SBATCH --ntasks=144                # must equal WRF_TASKS + FVCOM_TASKS below
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# #SBATCH --partition=<your_partition>
# #SBATCH --nodes=<N>               # set explicitly if you want WRF/FVCOM pinned to
#                                    # specific node counts rather than letting Slurm pack them

# ---------------------------------------------------------------------------
# Test/submit script for a coupled WRF-FVCOM-OASIS3-MCT run.
#
# It does three things, in order:
#   1. Pre-flight checks -- catch the failure modes that have come up while
#      setting this system up (mismatched task counts, drifted namcouple/grid
#      file copies, malformed namcouple, missing inputs) BEFORE burning
#      allocation time on a run that's going to abort in oasis_enddef anyway.
#   2. Launches wrf.exe and fvcom together as one MPMD job via
#      `srun --multi-prog`, each from its own working directory (each needs
#      its own namelist/tables vs. grid/sigma files, and each needs to see
#      namcouple + grids.nc/masks.nc locally).
#   3. Post-run checks -- did OASIS abort, did both sides actually produce
#      output.
#
# Edit the CONFIGURATION block below for your paths and task counts, then
# `sbatch submit_wrf_fvcom_coupled.sh`.
# ---------------------------------------------------------------------------

set -euo pipefail
source load_modules.sh

########################################
# 0. Configuration -- edit for your case
########################################
WRF_TASKS=120                        # must match the "awrf" count in namcouple $NBMODEL
FVCOM_TASKS=24                       # must match the "lfvc" count in namcouple $NBMODEL
TOTAL_TASKS=$((WRF_TASKS + FVCOM_TASKS))

WRF_RUNDIR=/compass/glm200001/cmu/coupled-run/nu-wrf-v11_cpl_oasis4/WRF/run
FVCOM_RUNDIR=/compass/glm200001/cmu/coupled-run/FVCOM41_oasis_wrf_fvcom_iceDynamic_new/run

WRF_EXE=wrf.exe
FVCOM_EXE=fvcom
FVCOM_CASE=gl                  # FVCOM case-name prefix, e.g. casename_grd.dat

LAUNCHER=mpirun

OMPI_MCA_coll_tuned_use_dynamic_rules=1
OMPI_MCA_coll_tuned_bcast_algorithm=1
OMPI_MCA_coll_hcoll_enable=0

########################################
# 1. Pre-flight checks
########################################
fail() { echo "PRE-FLIGHT FAIL: $1" >&2; exit 1; }
ok()   { echo "  OK: $1"; }

echo "== Pre-flight checks =="

# 1a. Requested Slurm allocation matches WRF_TASKS + FVCOM_TASKS.
#     (Falls back to TOTAL_TASKS when run outside Slurm, so this only bites
#     when it's actually meaningful.)
if [[ "${SLURM_NTASKS:-$TOTAL_TASKS}" -ne "$TOTAL_TASKS" ]]; then
    fail "SLURM_NTASKS (${SLURM_NTASKS:-unset}) != WRF_TASKS+FVCOM_TASKS ($TOTAL_TASKS). Fix #SBATCH --ntasks or the counts above."
fi
ok "task count: $WRF_TASKS (WRF) + $FVCOM_TASKS (FVCOM) = $TOTAL_TASKS"

# 1b. Executables exist and are executable.
[[ -x "$WRF_RUNDIR/$WRF_EXE" ]]     || fail "$WRF_RUNDIR/$WRF_EXE missing or not executable"
[[ -x "$FVCOM_RUNDIR/$FVCOM_EXE" ]] || fail "$FVCOM_RUNDIR/$FVCOM_EXE missing or not executable"
ok "wrf.exe and fvcom executables present"

# 1c. namcouple + grid description files present in BOTH run directories,
#     and NOT drifted copies of each other (should be symlinks to one master).
for f in namcouple grids.nc masks.nc; do
    [[ -e "$WRF_RUNDIR/$f" ]]   || fail "$f missing from $WRF_RUNDIR"
    [[ -e "$FVCOM_RUNDIR/$f" ]] || fail "$f missing from $FVCOM_RUNDIR"
    if ! cmp -s "$WRF_RUNDIR/$f" "$FVCOM_RUNDIR/$f"; then
        fail "$f differs between $WRF_RUNDIR and $FVCOM_RUNDIR -- symlink both to one master copy instead of maintaining separate copies"
    fi
done
ok "namcouple/grids.nc/masks.nc present and identical in both run dirs"

# 1d. namcouple $NBMODEL task counts match what we're actually launching.
#     ($NBMODEL's content is the line immediately following the keyword.)
nbmodel_line=$(grep -A1 '\$NBMODEL' "$WRF_RUNDIR/namcouple" | tail -1)
echo "  namcouple \$NBMODEL line: $nbmodel_line"
if ! echo "$nbmodel_line" | grep -qw "$WRF_TASKS" || ! echo "$nbmodel_line" | grep -qw "$FVCOM_TASKS"; then
    fail "namcouple \$NBMODEL task counts don't match WRF_TASKS=$WRF_TASKS / FVCOM_TASKS=$FVCOM_TASKS"
fi
ok "namcouple \$NBMODEL task counts match"

# 1e. $NFIELDS matches the number of field blocks actually present.
#     (Rough structural check, not a full namcouple parse -- catches the
#     "declared 11 fields but only wrote 9 blocks" class of error before
#     the OASIS parser does.)
nfields_declared=$(grep -A1 '\$NFIELDS' "$WRF_RUNDIR/namcouple" | tail -1 | tr -d '[:space:]')
nfields_found=$(grep -cE '(EXPORTED|EXPOUT|AUXILARY|IGNORED|IGNOUT|INPUT|OUTPUT)[[:space:]]*$' "$WRF_RUNDIR/namcouple" || true)
echo "  \$NFIELDS declared: $nfields_declared, field status lines found: $nfields_found"
if [[ "$nfields_declared" -ne "$nfields_found" ]]; then
    fail "\$NFIELDS ($nfields_declared) doesn't match number of field blocks found ($nfields_found) -- namcouple is likely malformed"
fi
ok "\$NFIELDS matches field block count"

# 1f. Required WRF inputs staged.
for f in wrfinput_d01 wrfbdy_d01 namelist.input; do
    [[ -e "$WRF_RUNDIR/$f" ]] || fail "$f missing from $WRF_RUNDIR"
done
ok "wrfinput_d01 / wrfbdy_d01 / namelist.input present"

# 1g. Required FVCOM inputs staged (adjust suffix list to your case).
for suffix in grd dep sigma obc cor; do
    [[ -e "$FVCOM_RUNDIR/input/${FVCOM_CASE}_${suffix}.dat" ]] || fail "${FVCOM_CASE}_${suffix}.dat missing from $FVCOM_RUNDIR"
done
[[ -e "$FVCOM_RUNDIR/${FVCOM_CASE}_run.nml" ]] || fail "${FVCOM_CASE}_run.nml missing from $FVCOM_RUNDIR"
ok "FVCOM grid/bathymetry/sigma/obc/namelist files present"

echo "== Pre-flight checks passed =="
echo

########################################
# 2. Build the MPMD launch config
########################################
MPMD_CONF=$(mktemp ./mpmd_conf.XXXXXX)
{
    echo "0-$((WRF_TASKS - 1)) $WRF_RUNDIR/wrf_wrapper.sh"
    echo "$WRF_TASKS-$((TOTAL_TASKS - 1)) $FVCOM_RUNDIR/fvcom_wrapper.sh"
} > "$MPMD_CONF"

# Wrapper scripts so each component runs from its own directory (its own
# namelist/tables vs. grid/sigma files) while still sharing one MPMD job.
cat > "$WRF_RUNDIR/wrf_wrapper.sh" <<EOF
#!/bin/bash
cd "$WRF_RUNDIR"
exec ./$WRF_EXE
EOF
chmod +x "$WRF_RUNDIR/wrf_wrapper.sh"

cat > "$FVCOM_RUNDIR/fvcom_wrapper.sh" <<EOF
#!/bin/bash
cd "$FVCOM_RUNDIR"
exec ./$FVCOM_EXE --CaseName $FVCOM_CASE --dbg_lvl=4 --dbg_par 
EOF
chmod +x "$FVCOM_RUNDIR/fvcom_wrapper.sh"

echo "== Generated MPMD config ($MPMD_CONF) =="
cat "$MPMD_CONF"
echo

########################################
# 3. Launch
########################################
echo "== Launching coupled run: $WRF_TASKS WRF tasks + $FVCOM_TASKS FVCOM tasks (launcher: $LAUNCHER) =="
set +e
if [[ "$LAUNCHER" == "srun" ]]; then
    srun --multi-prog "$MPMD_CONF"
elif [[ "$LAUNCHER" == "mpirun" ]]; then
    mpirun -np "$WRF_TASKS" "$WRF_RUNDIR/wrf_wrapper.sh" : \
           -np "$FVCOM_TASKS" "$FVCOM_RUNDIR/fvcom_wrapper.sh"
else
    echo "Unknown LAUNCHER='$LAUNCHER' (expected srun or mpirun)" >&2
    exit 1
fi
RUN_STATUS=$?
set -e

########################################
# 4. Post-run checks
########################################
echo
echo "== Post-run checks =="

if [[ $RUN_STATUS -ne 0 ]]; then
    echo "  WARNING: $LAUNCHER exited with status $RUN_STATUS"
fi

if grep -RIl "ABORT" "$WRF_RUNDIR"/debug.root.* "$FVCOM_RUNDIR"/debug.root.* 2>/dev/null; then
    echo "  WARNING: OASIS ABORT string found in the debug log(s) listed above"
else
    echo "  OK: no OASIS ABORT string found in debug logs"
fi

if ls "$WRF_RUNDIR"/wrfout_d01_* >/dev/null 2>&1; then
    echo "  OK: WRF output files present"
else
    echo "  WARNING: no wrfout_d01_* files found -- WRF likely didn't reach its first output step"
fi

if ls "$FVCOM_RUNDIR/${FVCOM_CASE}"_*.nc >/dev/null 2>&1; then
    echo "  OK: FVCOM output files present"
else
    echo "  WARNING: no FVCOM output NetCDF files found"
fi

echo
echo "== Done. For failures, check rsl.error.0000 (WRF), the FVCOM stdout/log, and debug.root.01 (OASIS) in each run dir. =="

