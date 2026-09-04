# WRF–FVCOM Coupling (Great Lakes)

Coupled **WRF** (atmosphere) ↔ **FVCOM** (ocean/lake, with ice dynamics) modelling system for
the Great Lakes, coupled through **OASIS3‑MCT**, plus the automation that keeps a long
multi‑year run advancing on the COMPASS HPC cluster: submitting restart cycles, downloading
ERA5 forcing, and rebuilding WRF boundary conditions with WPS.

---

## Current crontab schedule

```cron
0 8 * * * cd /compass/glm200001/cmu/coupled-run && sbatch scripts/glm_restart.sh
```

Every day at **08:00** (cluster local time) cron changes into the run root
`/compass/glm200001/cmu/coupled-run` and submits `scripts/glm_restart.sh` to SLURM.
`glm_restart.sh` is the automation entry point: each invocation either advances the coupled
run by one restart cycle, kicks off the WPS / ERA5 pipeline if boundary data is running out,
or notices there is nothing to do and exits cleanly. It self‑guards against overlapping
submissions, so a fixed daily trigger is safe.

---

## What's in this repository

| Path | Contents |
|------|----------|
| `nu-wrf-v11_cpl_oasis4/` | **WRF source tree** — NU‑WRF v11 with the OASIS3 coupling interface (`frame/module_cpl.F`, `frame/module_cpl_oasis3.F`). Also contains `WRF/run/` (run directory: `namelist.input`, `wrfrst_d01_*`, `wrfbdy_d01`, `wrfinput_d01`) and `WPS/` (ungrib / metgrid / real preprocessing, `namelist.wps`). |
| `FVCOM41_oasis_wrf_fvcom_iceDynamic_new/` | **FVCOM source tree** — FVCOM 4.1 with the OASIS coupling module `OASIS3MCT.F` (`module mod_var_cpl`) and the finite‑volume sea‑ice **dynamics** additions. Its `run/` directory holds `gl_run.nml`, `input/`, and `output/` (`gl_restart_*.nc`, history files). |
| `ERA5_download/` | **ERA5 forcing acquisition** via the Copernicus CDS API. `cdsapi-levels.py` (monthly pressure‑level requests), `cdsapi-surface.py` (annual surface requests), `submit_levels.sh` / `submit_surface.sh` (SLURM submit wrappers), `wget_levels.py` / `wget_surface.py` (collect finished downloads), `cdsapi_requests.csv` (request/download ledger). GRIB lands in `plevs-ERA5/plevs-ERA5-YYYY/` and `surface-ERA5/`. |
| `scripts/` | **Restart automation.** `glm_restart.sh` (cron driver), `check_wrf_inputs.sh` (boundary‑coverage check), `submit_WPS.sh` (rebuild WRF inputs), `submit_ERA5_download.sh` (fire CDS requests for a year), `wget_cdsapi_requests.sh` (collect ERA5 downloads), `README.md`. |
| `namcouple`, `namcouple.archive` | **OASIS3‑MCT configuration** — field list, coupling periods, lags, remapping (SCRIPR `DISTWGT` / `BILINEAR`), and restart‑file names for each exchanged field. `.archive` is a kept previous version. |
| `grids.nc`, `masks.nc` | **OASIS grid / mask description files** for the `fvcomN` (cell), `fvcomM` (node) and WRF (`wrf1_d01` / `wrf2_d01`) grids, written by the models on first run and reused thereafter. |
| `load_modules.sh` | HPC **environment**: purges and loads the compiler / MPI / HDF5 / NetCDF / Jasper stack and exports `NETCDF`, `HDF5`, `WRF_DIR`, `OASIS_DIR`, etc. Sourced by the scripts. |
| `submit_coupledrun.sh` | **SLURM launcher** for the coupled MPMD job (WRF + FVCOM + OASIS). Job name `coupled_run_glm`, account `glm200001`. |
| `small_submit_coupledrun.sh` | Reduced‑size / short‑test variant of the launcher. |

---

## The automated restart workflow (`scripts/glm_restart.sh`)

Each cron invocation runs one cycle:

1. **Guard against stacking** — if a coupled job (or a WPS job) is already queued/running, log and exit.
2. **Determine the restart time** — read the newest `wrfrst_d01_*` timestamp and the time in the
   newest `FVCOM_RUN/output/gl_restart_*.nc`; if they disagree, use the earlier of the two.
3. **Stage the FVCOM restart** — copy the chosen `output/gl_restart_*.nc` into `FVCOM_RUN/input/`
   so it can be referenced as `STARTUP_FILE`.
4. **Check boundary coverage** — `check_wrf_inputs.sh <restart_time>` reads `Times` from
   `wrfbdy_d01`. Exit 0 = still covers the next window; exit 10 = a WPS rebuild was launched.
5. **Branch:**
   - **Covered:** re‑arm `namelist.input` and `gl_run.nml` for a **hotstart** at the restart
     time, refresh the OASIS coupling restart files, and `sbatch submit_coupledrun.sh`.
   - **Not covered:** `submit_WPS.sh <restart_time>` builds fresh `wrfinput_d01` / `wrfbdy_d01`
     out to `WPS_WINDOW_MONTHS` (default 48), first verifying the needed ERA5 GRIB exists and
     firing `submit_ERA5_download.sh` for any missing years.
6. **Collect downloads** — `wget_cdsapi_requests.sh` pulls any finished ERA5 requests and
   updates `cdsapi_requests.csv`.
7. **Exit** — cron fires again on the next schedule.

### Hotstart hand‑off (per cycle)

- **FVCOM** (`gl_run.nml`): `STARTUP_TYPE = 'hotstart'`, `STARTUP_FILE` pointed at the staged
  `gl_restart_*.nc`, `STARTUP_{UV,TS,TURB}_TYPE = 'set values'`, `START_DATE` / `*_FIRST_OUT`
  set to the restart instant, `END_DATE` to the next cycle end.
- **WRF** (`namelist.input`): `restart = .true.`, `start_*` set to the restart time, run length
  to the cycle length.
- **OASIS**: `namcouple` unchanged; the timestamped coupling‑restart files
  (`TC<seconds>_rst_*.nc`) that match the restart instant are copied onto the base names
  (`rst_fvc_sst.nc`, `rst_wrf_airt.nc`, …) that OASIS reads at its time 0. FVCOM's coupling
  clock (`itap_sec`) restarts from 0 on a hotstart, so the restart instant must line up with a
  time for which both sides have a `TC*_rst_*.nc` file.

---

## `scripts/` reference

| Script | Args | Role |
|--------|------|------|
| `glm_restart.sh` | none | Cron driver. Job‑stacking guard, restart‑time detection, namelist re‑arm, submits the coupled run or the WPS pipeline, collects ERA5 downloads. Exports `GLM_COUPLED_ROOT`, writes `log.glm_restart`. |
| `check_wrf_inputs.sh` | `<YYYY‑MM‑DD_HH:MM:SS>` | Reads `Times` from `wrfbdy_d01` with `ncdump`; exit 0 if the boundary file still covers the upcoming window, exit 10 if it triggered `submit_WPS.sh`. |
| `submit_WPS.sh` | `<YYYY‑MM‑DD_HH:MM:SS>` | Rebuilds `wrfinput_d01` / `wrfbdy_d01` from the restart time out to `WPS_WINDOW_MONTHS` months. Confirms required ERA5 monthly pressure‑level and yearly surface GRIB exist (launches downloads if not), then runs ungrib → metgrid → real as dependent SLURM jobs and widens the namelist windows. |
| `submit_ERA5_download.sh` | `--input-year YYYY` | Fires asynchronous CDS requests for one calendar year — `cdsapi-levels.py` monthly, `cdsapi-surface.py` yearly — and appends rows to `cdsapi_requests.csv` with `pending` status. |
| `wget_cdsapi_requests.sh` | none | Walks `cdsapi_requests.csv` for pending rows, queries CDS job status, atomically downloads finished GRIB into `plevs-ERA5/…` / `surface-ERA5/`, marks rows complete, logs to `log.cdsapi_downloads`. Skips not‑ready jobs silently. |

---

## ERA5 forcing download

`cdsapi_requests.csv` is the ledger: one row per submitted CDS request with the year/month, a
CDS request id, and a `downloaded` flag. `submit_ERA5_download.sh` adds rows; `wget_cdsapi_requests.sh`
(and `ERA5_download/wget_levels.py` / `wget_surface.py`) close them out. Requests are split into
**monthly pressure‑level** GRIB (`era5_plevs_YYYY-MM.grib`) and **yearly surface** GRIB
(`era5_surf_YYYY.grib`) so a single failed month doesn't force a whole‑year re‑request.

Prerequisites: a working `~/.cdsapirc` (valid `url:` and `key:`) and the `cdsapi` Python
package (optionally in a dedicated venv — see environment variables).

---

## WPS processing

WPS is only run when `check_wrf_inputs.sh` finds `wrfbdy_d01` no longer reaches the next
coupling window. `submit_WPS.sh` then:

1. Computes the new window `[restart_time, restart_time + WPS_WINDOW_MONTHS]`.
2. Checks the ERA5 GRIB inventory for that span; triggers `submit_ERA5_download.sh` for missing years and stops (a later cron cycle resumes once GRIB has arrived).
3. Runs **ungrib → metgrid → real** as chained SLURM jobs against `namelist.wps` / `namelist.input`.
4. Produces fresh `wrfinput_d01` and `wrfbdy_d01` and widens the WRF/FVCOM run windows.

---

## Coupling configuration

- **`namcouple`** — the OASIS3‑MCT field table. Ocean→atmosphere sends SST and ice
  concentration (`fvcom_send_sst`, `fvcom_send_aice`) from the FVCOM node grid (`fvcomM`) to
  WRF (`DISTWGT`); atmosphere→ocean sends 2 m air temperature, radiation, humidity, winds,
  pressure, precip/evap (`WRF_d01_EXT_d01_*` → `fvcom_recv_*`) to the FVCOM cell grid
  (`fvcomN`, `BILINEAR`). Coupling period and lag are the `dt` / `LAG` columns; the 6th column
  is the restart‑file base name.
- **`grids.nc` / `masks.nc`** — grid coordinates and land/sea masks for every OASIS grid,
  written by the models on the first run.
- **Coupling restart cadence** — WRF (`module_cpl_oasis3.F`, `cpl_oasis_snd`) and FVCOM
  (`OASIS3MCT.F`, `tocoupler`) each decide when `oasis_put` writes a `TC<seconds>_rst_*.nc`
  file. These files are what the restart automation copies onto the base names for the next
  hotstart. (See **Open questions** for the interval choice.)
- **FVCOM never calls `oasis_terminate`** — it shuts down through `MPI_FINALIZE`, so the only
  FVCOM‑side coupling restarts are the mid‑run `TC*_rst_fvc_*.nc` files.

---

## Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `GLM_COUPLED_ROOT` | parent of `scripts/` | Project root override; exported by `glm_restart.sh`. |
| `GLM_LOG` | `$ROOT/log.glm_restart` | Main timestamped log. |
| `GLM_PYTHON` | `python3` | Interpreter for the CDS API scripts. |
| `GLM_CDSAPI_VENV` | unset | Optional path to a venv containing `cdsapi`. |
| `CDSAPI_RC` | `$HOME/.cdsapirc` | CDS credentials file (parsed for `url:` / `key:`). |
| `WPS_WINDOW_MONTHS` | `48` | Length of each rebuilt boundary window. |

Path resolution in every script:

```bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${GLM_COUPLED_ROOT:-$(dirname "$SCRIPT_DIR")}"
```

Key directories (relative to `ROOT`):

```
ROOT/
├── scripts/
├── load_modules.sh
├── namcouple                 grids.nc  masks.nc
├── submit_coupledrun.sh
├── nu-wrf-v11_cpl_oasis4/
│   ├── WRF/run/              $WRF_RUN   (wrfrst_d01_*, wrfbdy_d01, namelist.input)
│   └── WPS/                  $WPS_DIR   (ungrib/metgrid, namelist.wps)
├── FVCOM41_oasis_wrf_fvcom_iceDynamic_new/
│   └── run/                  $FVCOM_RUN (input/, output/gl_restart_*.nc, gl_run.nml)
└── ERA5_download/
    ├── plevs-ERA5/plevs-ERA5-YYYY/    (era5_plevs_YYYY-MM.grib)
    ├── surface-ERA5/                  (era5_surf_YYYY.grib)
    └── cdsapi_requests.csv
```

---

## Running by hand

```bash
cd /compass/glm200001/cmu/coupled-run
scripts/glm_restart.sh            # one cycle, same as cron
# or force a coupled run of the current namelists:
sbatch submit_coupledrun.sh
```

Monitor:

```bash
tail -f /compass/glm200001/cmu/coupled-run/log.glm_restart
tail -f /compass/glm200001/cmu/coupled-run/log.cdsapi_downloads
tail -f /compass/glm200001/cmu/coupled-run/nu-wrf-v11_cpl_oasis4/WRF/run/rsl.out.0000
squeue -A glm200001
```

