# `scripts/` — coupled WRF/FVCOM Great Lakes run automation

These scripts drive a long, self-restarting coupled **WRF + FVCOM** (OASIS3-MCT)
simulation of the Great Lakes. The run is executed as a chain of short SLURM jobs:
each cycle restarts the model from the last set of restart files, and — when the
model catches up to the end of the prepared boundary window — a new block of WRF
inputs is built from ERA5 first, then the run continues. The block length is set
by `WPS_WINDOW_MONTHS` at the top of `submit_WPS.sh` (default 48 months).

`glm_restart.sh` is the single entry point. It is designed to be run **repeatedly
from cron**; every invocation either advances the run by one cycle, or notices
there is nothing to do and exits cleanly.

```
cron ──► glm_restart.sh ──┬─► check_wrf_inputs.sh ──► submit_WPS.sh ──► (sbatch) submit_metgrid.sh
                          │                                        └──► (sbatch) submit_real.sh
                          │                                        └──► submit_ERA5_download.sh ──► cdsapi-levels.py / cdsapi-surface.py
                          ├─► (sbatch) submit_coupledrun.sh
                          └─► wget_cdsapi_requests.sh
```

All of these scripts live in this directory except `load_modules.sh`, which stays
in the coupled-run root (it is also sourced by cluster job scripts by absolute
path).

---

## Layout / path convention

Every script resolves its own location and the project root the same way:

```sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../coupled-run/scripts
ROOT="${GLM_COUPLED_ROOT:-$(dirname "$SCRIPT_DIR")}"          # .../coupled-run
```

* Sibling scripts are called via `"$SCRIPT_DIR/..."`.
* Everything else (model run dirs, `ERA5_download/`, logs, `load_modules.sh`) is
  under `"$ROOT/..."`.
* `glm_restart.sh` `export`s `GLM_COUPLED_ROOT`, so child scripts inherit a
  consistent root. Set `GLM_COUPLED_ROOT` yourself to run any script against a
  different checkout.

Referenced model directories:

| Variable    | Path                                             |
|-------------|--------------------------------------------------|
| `WRF_RUN`   | `$ROOT/nu-wrf-v11_cpl_oasis4/WRF/run`            |
| `WPS_DIR`   | `$ROOT/nu-wrf-v11_cpl_oasis4/WPS`                |
| `FVCOM_RUN` | `$ROOT/FVCOM41_oasis_wrf_fvcom_iceDynamic_new/run` |
| `ERA5_DIR`  | `$ROOT/ERA5_download`                            |

---

## Environment variables

| Variable            | Default                     | Used by | Purpose |
|---------------------|-----------------------------|---------|---------|
| `GLM_COUPLED_ROOT`  | parent of `scripts/`        | all     | Project root override. Exported by `glm_restart.sh`. |
| `GLM_LOG`           | `$ROOT/log.glm_restart`     | all     | Main log file. Exported by `glm_restart.sh`. |
| `GLM_PYTHON`        | `python3`                   | `submit_ERA5_download.sh` | Python interpreter for the cdsapi scripts. |
| `GLM_CDSAPI_VENV`   | *(unset)*                   | `submit_ERA5_download.sh` | If set and `$GLM_CDSAPI_VENV/bin/activate` exists, that venv is activated before calling the cdsapi scripts. |
| `CDSAPI_RC`         | `$HOME/.cdsapirc`           | `wget_cdsapi_requests.sh` | CDS API credentials file (parsed for `url:` and `key:`). |

---

## Logs (written in `$ROOT`)

| File                    | Written by | Contents |
|-------------------------|------------|----------|
| `log.glm_restart`       | every script here | Timestamped, human-readable trace of each cycle. Each line is prefixed `[YYYY-MM-DD HH:MM:SS]` and, for sub-scripts, a tag (`check_wrf_inputs:`, `submit_WPS:`, `submit_ERA5_download:`). |
| `log.cdsapi_downloads`  | `wget_cdsapi_requests.sh` | One line per successfully downloaded ERA5 grib. |

---

## Scripts

### `glm_restart.sh` — top-level driver

**Usage:** `./glm_restart.sh` (no arguments). Intended to run from cron on a
fixed interval.

Each invocation:

1. **Don't stack jobs.** Runs `squeue -u $USER` and greps job names
   `coupled_run_glm`, `real_glm`, `metgrid_glm`. If any is queued/running, it
   logs that and exits 0 immediately — no further action.
2. **Find the restart point.**
   * WRF: newest `wrfrst_d01_*` in `WRF_RUN`; timestamp parsed from the filename.
   * FVCOM: newest `output/gl_restart_*.nc` in `FVCOM_RUN`; timestamp read from
     the last `Times` entry via `ncdump`.
   * If the two disagree, it logs the discrepancy and proceeds with the
     **earliest** of the two.
3. **Validate WRF inputs** by calling `check_wrf_inputs.sh <restart_ts>`:
   * **exit 0** — boundary data still covers the upcoming window. Re-arm the
     namelists for a restart and submit the coupled run:
     - `WRF_RUN/namelist.input`: set `start_*` to the restart time, `end_*` to
       restart + run length, `run_hours/minutes/seconds` to that length,
       `run_days = 0`, `restart = .true.`. Run length is read from the current
       `run_days/hours/minutes/seconds` in the namelist (fallback: 43200 s = 12 h,
       the namcouple `$RUNTIME`).
     - `FVCOM_RUN/gl_run.nml`: set `START_DATE`/`END_DATE`, `STARTUP_TYPE =
       'hotstart'`, `STARTUP_FILE` to the latest `gl_restart_*.nc`, and the
       `RST/NC/NCAV_FIRST_OUT` times to the restart time.
     - `sbatch --parsable submit_coupledrun.sh` (job name `coupled_run_glm`).
   * **exit 10** — `check_wrf_inputs.sh` launched the WPS pipeline. Nothing else
     to do this cycle; the coupled run resumes on a later cycle once `real_glm`
     finishes.
   * **any other code** — logged as unexpected; no coupled run submitted.
4. **Always** runs `wget_cdsapi_requests.sh` at the end, regardless of the branch
   taken above (via the `finish()` helper, which also logs the end banner).

**Reads:** `wrfrst_d01_*`, `gl_restart_*.nc`, `namelist.input`, `gl_run.nml`.
**Writes:** `namelist.input`, `gl_run.nml`, `log.glm_restart`. **Submits:**
`coupled_run_glm`.

---

### `check_wrf_inputs.sh` — is the boundary window still valid?

**Usage:** `./check_wrf_inputs.sh <YYYY-MM-DD_HH:MM:SS>`

Reads the `Times` variable from `WRF_RUN/wrfbdy_d01` with `ncdump`. If the
restart timestamp is **at or after the second-to-last** boundary time, the
existing `wrfbdy_d01` no longer covers the coming coupling window, so it calls
`submit_WPS.sh <restart_ts>` to rebuild inputs.

`wrfbdy_d01` missing, or fewer than 2 `Times` entries, also triggers WPS.

If `ncdump` is not on `PATH`, it first tries `source "$ROOT/load_modules.sh"`.

**Exit codes:**

| Code | Meaning |
|------|---------|
| `0`  | Inputs still valid — no WPS needed. |
| `10` | WPS pipeline was launched (via `submit_WPS.sh`). |
| `2`  | Usage error or an unparseable timestamp. |

---

### `submit_WPS.sh` — rebuild WRF inputs for a new boundary window

**Usage:** `./submit_WPS.sh <YYYY-MM-DD_HH:MM:SS>`

Builds `wrfinput_d01` / `wrfbdy_d01` for a fresh window running from
`<restart_ts>` (rounded down to the hour) to **restart + `WPS_WINDOW_MONTHS`**.

`WPS_WINDOW_MONTHS` is a variable at the top of the script (default `48`, i.e.
4 years). Set it to change the length of the window that gets built; the value is
in **months** and every downstream date (namelist windows, ERA5 year/month
coverage, `run_days`) is derived from it.

1. **ERA5 availability.** Checks that every
   `ERA5_download/plevs-ERA5/plevs-ERA5-YYYY/era5_plevs_YYYY-MM.grib` (each month
   in the window) and `ERA5_download/surface-ERA5/era5_surf_YYYY.grib` (each year)
   exists. If anything is missing, it calls `submit_ERA5_download.sh
   --input-year YYYY` for every year in the window and **aborts this cycle**
   (exit 1) — downloads are asynchronous, so a later cycle picks up where this
   one left off.
2. `cd "$WPS_DIR"`, `source "$ROOT/load_modules.sh"`.
3. Sets `start_date` / `end_date` in `namelist.wps`.
4. **Ungrib pressure levels:** `prefix = 'ERA5-PRES'`, `link_grib.csh` the
   monthly plevs gribs, run `ungrib.exe`, verify `ungrib.log` reports
   "Successful completion of ungrib" and that `ERA5-PRES:*` intermediates exist.
5. **Ungrib single levels:** same, with `prefix = 'ERA5-SURF'` over the yearly
   surface gribs.
6. `sbatch --parsable submit_metgrid.sh <start> <end>` → job `metgrid_glm`
   (produces `met_em.*` for the window).
7. `sbatch --parsable --dependency=afterok:<metgrid> submit_real.sh <start> <end>`
   → job `real_glm` (refreshes the `met_em` links in `WRF/run`, runs `real.exe`,
   produces `wrfinput_d01` / `wrfbdy_d01`).
8. **Widen the namelists** to the new window: `WRF/run/namelist.input`
   `start_*`/`end_*`, `run_days` = window length in days,
   `run_hours/minutes/seconds = 0`, `restart = .false.`; and
   `FVCOM_RUN/gl_run.nml` `END_DATE` to the new end time.

Every step is verified before the next runs; failures call `die` (log + exit 1).
Returns 0 once the batch jobs are queued.

**Note:** while `metgrid_glm` / `real_glm` are queued, `glm_restart.sh` stops at
its `squeue` guard, so nothing overlaps. Once `real_glm` finishes and
`wrfbdy_d01` covers the window, the next `glm_restart.sh` submits the coupled run.

---

### `submit_ERA5_download.sh` — request one year of ERA5 forcing

**Usage:** `./submit_ERA5_download.sh --input-year YYYY`

Fires **asynchronous** CDS requests for one calendar year:

* pressure levels — one request per month (up to 12), via `cdsapi-levels.py`;
* single levels — one request for the whole year, via `cdsapi-surface.py`.

The Python scripts run with `wait_until_complete=False` and print
`REQUEST_ID=<id>` (or `EXISTS` / `ERROR`); this script parses that line. Nothing
is downloaded here — `wget_cdsapi_requests.sh` collects finished files on a later
cycle.

A request is **skipped** when the target grib already exists on disk, or an
undownloaded row for the same year/month is already present in the CSV.

**Writes:** `ERA5_download/cdsapi_requests.csv` — appended, header created if
absent:

```
request_timestamp,request_year,request_month,request_id,downloaded
```

`request_month` is blank for surface (whole-year) requests; `downloaded` is `0`
(pending) or `1` (done, set later by `wget_cdsapi_requests.sh`).

**Exit codes:** `0` normal, `2` bad/missing `--input-year`.

If `GLM_CDSAPI_VENV` points at a venv with `bin/activate`, it is sourced first;
otherwise `GLM_PYTHON` (default `python3`) must provide the `cdsapi` package.

---

### `wget_cdsapi_requests.sh` — collect finished ERA5 downloads

**Usage:** `./wget_cdsapi_requests.sh` (no arguments). Run at the end of every
`glm_restart.sh` cycle.

Walks `ERA5_download/cdsapi_requests.csv`. For each row with `downloaded == 0`
(and a real `request_id` — not `ERROR` / `EXISTS`):

1. Derives the destination path from the row (`plevs-ERA5/...` if a month is set,
   else `surface-ERA5/...`).
2. If the file already exists, flips the flag to `1`.
3. Otherwise queries `"$CDS_URL/retrieve/v1/jobs/<id>/results"` (URL and key read
   from `CDSAPI_RC`, default `~/.cdsapirc`). If the job is finished, `wget`s the
   asset to `<dest>.part`, then atomically `mv`s it into place and sets
   `downloaded = 1`.
4. Jobs that are not ready / not available are **skipped silently** — no error,
   no warning; the row is left untouched for a future cycle.

The CSV is rewritten only if something changed. Successful downloads are logged
to `log.cdsapi_downloads`.

Exits 0 (a no-op) if the CSV or the credentials file is missing.

---

### `submit_coupledrun.sh` — the coupled-model SLURM job

The batch script `glm_restart.sh` submits with `sbatch` (job name
`coupled_run_glm`, account `glm200001`). It is launched by the driver, not run
directly, and reads nothing from this directory — it operates in the WRF/FVCOM
run dirs configured inside it.

> ⚠️ This file is currently a **stub** (SBATCH headers only). Fill in the
> `srun`/MPMD launch for the OASIS-coupled WRF + FVCOM executables before using
> the pipeline for a real run.

---

## Related files outside this directory

| File | Role |
|------|------|
| `../load_modules.sh` | `module purge` + loads the compiler/MPI/HDF5/NetCDF stack and exports `NETCDF`, `HDF5`, `JASPER*`, `WRF_DIR`, `OASIS_DIR`, etc. Sourced by `glm_restart.sh`, `check_wrf_inputs.sh`, `submit_WPS.sh`, and by the cluster `submit_real.sh` / `submit_metgrid.sh` (the latter two by absolute path). |
| `../ERA5_download/cdsapi-levels.py` | Submits one async CDS request for one month of pressure-level ERA5. Prints `REQUEST_ID=` / `TARGET=`. |
| `../ERA5_download/cdsapi-surface.py` | Submits one async CDS request for a full year of single-level ERA5. |
| `../ERA5_download/cdsapi_requests.csv` | Ledger of submitted requests and their download state (see above). |
| `../nu-wrf-v11_cpl_oasis4/WPS/submit_metgrid.sh` | `metgrid_glm` job — runs `metgrid.exe` for the window. `sbatch`ed by `submit_WPS.sh`. |
| `../nu-wrf-v11_cpl_oasis4/WRF/run/submit_real.sh` | `real_glm` job — refreshes `met_em` links and runs `real.exe`. `sbatch`ed by `submit_WPS.sh` with `--dependency=afterok` on metgrid. |

---

## Typical cron setup

```cron
# every 30 min: advance the coupled Great Lakes run
*/30 * * * *  /compass/glm200001/cmu/coupled-run/scripts/glm_restart.sh >/dev/null 2>&1
```

`glm_restart.sh` is idempotent per cycle and self-guards against overlapping
jobs, so a short interval is safe. Watch progress with:

```sh
tail -f /compass/glm200001/cmu/coupled-run/log.glm_restart
```
