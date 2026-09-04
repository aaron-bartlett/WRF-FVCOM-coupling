* **Atmosphere:** `nu-wrf-v11_cpl_oasis4/` — NU‑WRF v11, ARW core, OASIS3 coupling under CPP key `key_cpp_oasis3`.
* **Ocean/ice:** `FVCOM41_oasis_wrf_fvcom_iceDynamic_new/` — FVCOM 4.1 with the dynamic‑ice (CICE‑derived) option, OASIS coupling under CPP key `oasis_coupler`.
* **Coupler:** OASIS3‑MCT, driven by `namcouple` (`$NBMODEL  2  wrfexe  ocean1`).
---

## 1. How the coupling is wired (context for the changes below)

### 1.1 Exchanged fields (`namcouple`)

| Direction | Fields | OASIS grids | Interp |
|---|---|---|---|
| FVCOM → WRF | `SST`, `AICE` (ice fraction) | `fvcomM` (node, `MGL`) → `wrf2_d01` | DISTWGT |
| WRF → FVCOM | `airt`, `airp`, `rh`, `spq`, `cloud`, `u10`, `v10`, `shortwave`, `longwave` (and, when enabled, `evap`, `precip`, `latent`, `sensible`, `netheat`, `dew`) | `wrf1_d01` → `fvcomN` (element, `NGL`) | BILINEAR |
| FVCOM ↔ SWAN | `u`, `v`, `zeta` out; `fx`, `fy` (radiation‑stress gradient) or `hs/wlen/dir/rtp/tmbot/ubot` (wave params) in | `fvcomN` / `fvcomM` | — |

### 1.2 Runtime switches

* **WRF:** `coupler_on` (compile) + `coupler_name=='oasis'`. Fields default to `'not defined'`
  and are activated by their presence in `namcouple`.
* **FVCOM:** `USE_OASIS_COUPLER` (from `&NML_OASIS_COUPLER` in `<case>_run.nml`) plus per‑path
  sub‑switches `OASIS_atm2ocn`, `OASIS_ocn2atm`, `OASIS_ocn2wav`, `OASIS_wav2ocn_RadStr`,
  `OASIS_wav2ocn_Param`, and the ice‑control knobs `OASIS_ICE_CTRL_AICE` /
  `OASIS_ICE_CTRL_SST` / `OASIS_CTRL_VAL_*` / `OASIS_freezing_temp`.

### 1.3 The FVCOM coupling data path (recurring pattern)

Almost every FVCOM change below is one link in this chain:

```
OASIS oasis_get ─▶ field_fvcom_recv_<v>(:,1)        (OASIS3MCT.F : fromcoupler)
                └▶ <V>_cell2  (this coupling slice) ─┐
                   <V>_cell0  (previous slice)       ├─ time-interp in mod_force.F:
                                                     │    <V>_cell = NEXT_WGHT_oasis*<V>_cell2
                   FTM%NEXT_WGHT_oasis /             │             + PREV_WGHT_oasis*<V>_cell0
                       PREV_WGHT_oasis  ────────────┘
                └▶ E2N2D ─▶ T_AIR / RH_AIR / PA_AIR / DSW_AIR / … (native FVCOM forcing arrays)
                └▶ AEXCHANGE (halo fill across MPI partition boundaries)
```

The `fvcom.F` main loop does the OASIS `oasis_put`/`oasis_get` and fills `*_cell0/2`; `mod_force.F`
does the temporal interpolation and substitutes the result for the file‑based forcing;
`mod_nctools.F` maintains the interpolation weights.

---

## 2. WRF side — `nu-wrf-v11_cpl_oasis4` (`chenfu…` tags)

### 2.1 `frame/module_cpl.F` — coupling field catalogue & fill

| Line(s) | Tag | Change | Role in coupling |
|---|---|---|---|
| 30 | `chenfu2017` | `INTEGER,PUBLIC,save :: nsecrun_save` | Public copy of the coupling clock (seconds since job start). |
| 230 | `chenfu2017` | `nsecrun_save = nsecrun` inside `cpl_settime` | Lets the physics layer (surface driver) read the OASIS time index it must pass to `oasis_get`. |
| 79, 90 | `chenfu2017` | `rcvname(:,:,4) = …'AICE'`; `sndname(:,:,7..28)` = `evap, precip, latent, sensible, shortwave, longwave, netheat, spq, rh, u10, v10, airt, airp, cloud, ru10, rv10, ust, u10tau, v10tau, ru10tau, rv10tau, dew` | Extends the stock OASIS‑WRF field list (SST/UOCE/VOCE + flux bundle) with the full near‑surface met state FVCOM needs to run its own bulk‑flux/ice thermodynamics, plus **AICE received back** from FVCOM. |
| 192 | `chenfu2019` | commented `grid%cplmask = grid%landmask*(-1)+1` | Note that the coupling mask is the sea mask (1 on water); actual masking now done in `module_cpl_oasis3.F`. |
| 348–351 | `chenfu2019` | extra locals `uu1,vv1,wd1,wtau,uu2,vv2,wd2,icnt2` in `cpl_snd` | Scratch for the relative‑wind / wind‑stress send fields. |
| 414–415 | `chenfu2017` | `cpl_snd` block "added for wrf‑fvcom‑swan" | Fills the new send fields from `grid%t2-273.15` (→ °C), `grid%psfc`, `grid%q2`, `grid%cldt`, `grid%u10`, `grid%v10`, … |
| 527 | `chenfu2019` | `ifldid = cpl_get_fldid('ru10')`; `cplsnd = grid%u10 - grid%uoce` (and `rv10`) | **Relative** 10 m wind (air minus ocean‑surface current) — the physically correct forcing for a current‑aware ocean/wave model. |
| 722 | `chenfu202005` | `llmust_store = llmust_store .OR. (cpl_get_fldid('SST') .gt. 0)` (was `cpl_toreceive(...)`) | Store WRF's read‑in SST into `SST_INPUT` whenever an SST coupling field is defined, so the coupler‑supplied SST can overwrite `grid%sst` without losing the background field. |

### 2.2 `frame/module_cpl_oasis3.F` — OASIS partition / grid / put

| Line(s) | Tag | Change | Role |
|---|---|---|---|
| 120–122 | `chenfu2017` / `chenfu2019` | `il_flag`, `oasis_grid_name`, `id_partA(max_domains)`, `id_partB(max_domains)` | Per‑nest partition ids; a **second** partition id so send and receive can carry different masks. |
| 172 | `chenfu2019` | `oasis_def_partition(id_partA(pgrid%id), …)` | Box partition keyed by WRF nest id (multi‑domain safe). |
| 174–189 | `chenfu2017` | 2nd `oasis_def_partition(id_partB…)`; `oasis_start_grids_writing`/`oasis_write_grid`/`oasis_write_mask` for `wrf1_d<NN>` (WRF as **source**) and `wrf2_d<NN>` (WRF as **destination**), masks from `grid%cplmask` (sent grid uses `mask*0` = all active; received grid uses `1-mask`) | Writes `grids.nc`/`masks.nc` for SCRIP weight generation. The dual grid is why `namcouple` names `wrf1_d01`/`wrf2_d01` separately. |
| 219, 248 | `chenfu2019` | `oasis_def_var(… trim(name), id_partA(pgrid%id) …)` for every send/recv field | Trims blank‑padded names and binds each var to the nest partition. |
| 351–377 | `chenfu2019` | in `cpl_snd`: scan `namcouple` tables (`namsrcfld(iwd)(1:3)=='WRF'`) for `cpl_dt`/`cpl_lag`, set `irestart = mod(ksec+cpl_lag,cpl_dt)==0`, pass `write_restart=irestart` to `oasis_put` | WRF writes its OASIS restart (coupling‑field snapshot) exactly on coupling steps — needed for a clean coupled restart. Mirrors FVCOM's `tocoupler`. |

### 2.3 Other WRF files

| File : line | Tag | Change | Role |
|---|---|---|---|
| `frame/module_driver_constants.F:92` | `chenfu201902` | `max_cplfld = 30` (was 20) | Array bound for the enlarged send list (28 fields). |
| `phys/module_surface_driver.F:37–41, 875–879` | `chenfu2018` | thread `AICE` (`#ifdef key_cpp_oasis3`) through `surface_driver` dummy args + declaration `REAL … AICE INTENT(INOUT)` | Makes the FVCOM ice fraction visible inside the surface layer. |
| `phys/module_surface_driver.F:405` | `chenfu2017` | `USE module_cpl, ONLY : nsecrun_save` | Time index for the `cpl_rcv` calls. |
| `phys/module_surface_driver.F:1507` | `chenfu2017` | `CALL cpl_rcv(id,'AICE', … ,AICE)` alongside SST/UOCE/VOCE | Pulls ice fraction from OASIS each coupling step. |
| `phys/module_surface_driver.F:1518` | `chenfu201902` | `XICE = AICE` when `fractional_seaice == 1` | **Drives WRF's prognostic sea‑ice fraction `XICE` from FVCOM's ice model** — the core FVCOM→WRF feedback for ice. |
| `dyn_em/module_first_rk_step_part1.F:572–576` | `chenfu2018` | pass `AICE=grid%aice` into the `surface_driver` call (`#ifdef key_cpp_oasis3`) | Connects the Registry array to the surface driver. |
| `Registry/Registry.EM_COMMON:877` | `chenfu2017` | `state real AICE ij misc … "SEA SURFACE ICE COVER"` | Allocates/streams the received ice‑fraction field. |
| `arch/configure*.defaults`, `arch/postamble*`, `nu-wrf.cfg` | `--chenfu` (undated) | add `-Dkey_cpp_oasis3`, OASIS include/module dirs and link libs to the build | Build plumbing. *(Marked `--chenfu` without a year — listed for completeness.)* |

---

## 3. FVCOM side — `FVCOM41_oasis_wrf_fvcom_iceDynamic_new` (`ch…` tags)

There are ~320 tagged lines. They fall into the chronological families below. File:line
references are representative, not exhaustive (repeated identical edits — e.g. the ~30
`IF(ICE_MODEL) THEN !--ch2017 added` guards and the ~50 `!-ch2018_hrrr` reader alternates —
are described once with a count).

### 3.1 `ch2015` / `ch201508` — coupling clock & error‑check scaffolding

| File : line | What | Role |
|---|---|---|
| `internal_step.F:125` | `CHARACTER(LEN=80) :: SIM_TIME ! px201502 for coupling(ch2015)` | Wall‑clock string used to gate month‑based logic (deep‑water cap, monthly water level, monthly restart). |
| `internal_step.F:126` | `INTEGER :: IERR !ch201508` | MPI status for the added barriers/collectives. |
| `mod_main.F:934` | `CHARACTER(LEN=80) :: TIME_SAVE !--ch201508` | Saved model time string. |
| `mod_main.F:1116–1121` | `!----ch201508_error_check`: `CHECK_ERROR`, `CHECK_index(:)`, `CHECK_index_global(:)` | Global "something went non‑physical" flag, collected across ranks. Consumed by `nan_check` (3.9) and OASIS de‑overlap logic. |
| `mod_ncdio.F:199` | `!--ch201508 added for output data before error/nan` | Force an output dump just before a detected blow‑up, for debugging coupled runs. |

### 3.2 `ch201610` / `ch201610_oasis` — OASIS bootstrap

| File : line | What | Role |
|---|---|---|
| `mod_main.F:111` (`CONTROL`) | `il_commlocal`, `id_comp_fvcom`, `LOGICAL :: USE_OASIS_COUPLER` | The MPI communicator OASIS hands back to FVCOM (a split of `MPI_COMM_WORLD`), the OASIS component id, and the master runtime switch. |
| `mod_utils.F:274, 303` | `use CONTROL, only: il_commlocal` | Utility barriers/broadcasts must use the **coupled** communicator, not `MPI_COMM_WORLD` (which now also contains WRF ranks). |
| `mod_par.F:35`, `mod_par_special.F:35,250,260` | `!-----ch201610_oasis` … `IF (USE_OASIS_COUPLER)` | Decomposition/`flag_coupler` handling in the parallel setup. |
| `fvcom.F:507` | `use mod_var_cpl / mod_oasis` | Pull in the coupler module (`OASIS3MCT.F : module mod_var_cpl`). |
| `fvcom.F:756–856` | `ch201610 added for coupler` — `call ini_coupler('ocean1')` + cold‑start initial exchange | One‑time OASIS init: partitions, grid write, `oasis_def_var` for every `fvcom_send_*`/`fvcom_recv_*`, `oasis_enddef`, then the **first** send/receive so neither model deadlocks at `t=0` (see also `ba2026` reorder note, §4). |
| `fvcom.F:947–1345` | `ch201610 added for coupler` — the per‑timestep coupler block wrapping `INTERNAL_STEP` | Computes `itap_sec`, averages SST/AICE over the period, `tocoupler`/`fromcoupler` for atm and wave paths, deposits received data into `*_cell0/2`, halo‑exchanges. |
| `OASIS3MCT.F_ori:1109–1124` | `ch201610 added for corner dimension` | Grid‑corner arrays for `CONSERV` interpolation. *(In `OASIS3MCT.F_ori`, a retained backup of the coupler module; the live build uses `OASIS3MCT.F`, and `namcouple` notes CONSERV was abandoned because corners are never actually written.)* |
| `mod_petsc.F:39, 235` | `use control, only: il_commlocal`; `PETSC_COMM_WORLD = il_commlocal` | The semi‑implicit solver must run on the FVCOM‑only sub‑communicator. |

### 3.3 `ch201611` / `ch201611_oasis` — the bulk of the FVCOM implementation

**Dual OASIS partition/grid.** `OASIS3MCT.F` defines two partitions: `part_id` over elements
(`fvcomN`, `NGL` — all atm→ocean fields land here) and `part_id2` over nodes (`fvcomM`,
`MGL` — SST/AICE out, wave params in). Nodes shared across MPI partitions are de‑overlapped
(`local_nn`, `check_index_global`) so each global node is owned once.

**Received‑field buffering + time interpolation.**

| File : line | What | Role |
|---|---|---|
| `mod_ncll.F:136–137` | add `PREV_WGHT_oasis`, `NEXT_WGHT_oasis` to the `NCFTIME` type | Storage for the coupler‑specific temporal weights, parallel to the file weights `PREV_WGHT/NEXT_WGHT`. |
| `mod_nctools.F:9439, 9606` | `denom_oasis/numer_oasis`; compute `NEXT_WGHT_oasis = mod(seconds(NOW),cpl_dt_atm2ocn)/cpl_dt_atm2ocn` in `UPDATE_FILE_BRACKET` | Linear‑in‑time weight between the previous and next coupling slice. *(Relocated to the top of the routine by `ba2026` — §4.)* |
| `mod_main.F:1370–2172` | `ch201611_oasis` — `UUWIND_oasis{,0,2}`, `VVWIND_oasis{,0,2}`, `UUSTRX_oasis*`, `VVSTRX_oasis*`, `QPREC_cell{,0,2}`, `QEVAP_cell{,0,2}`, `QPREC2_cell*`, `QEVAP2_cell*`, `sst_on_cell/node`, `aice_on_cell/node`, `sst_sum/avg`, `aice_sum/avg`, `oasis_check` + their `ALLOCATE`s | The `_cell0` (previous), `_cell2` (next), `_cell` (interpolated) triplets for every coupled variable; SST/AICE period accumulators. |
| `mod_force.F:12472, 12514, 12743, 12892, 13414, 13472, 14013, 14056, 14344, 14388, 14984, 15073` | `ch201611 added for oasis` — in `UPDATE_HEAT_*` and `UPDATE_WIND`: `IF (USE_OASIS_COUPLER .and. OASIS_atm2ocn)` then `T_AIR_cell = NEXT_WGHT_oasis*T_AIR_cell2 + PREV_WGHT_oasis*T_AIR_cell0` (likewise RH/DLW/DSW/PA/HEAT_NET/UUWIND/VVWIND), `E2N2D` to node arrays, `UUWIND=UUWIND_oasis` | **This is where OASIS data replaces the file forcing.** The rest of FVCOM's heat/wind/pressure physics is untouched — it just sees `T_AIR`, `DSW_AIR`, `WSTRX`, … as if read from a file. |
| `mod_force.F:489, 14861–15105` | `ch201611 added if defined for ice` | Same override, `ICE_FILE` branch, so the dynamic‑ice thermodynamics is coupler‑driven too. |
| `mod_ncdio.F:3689–3812`, `mod_nctools.F` | `ch201611 added for oasis` — allocate + write `u10_oasis` / `v10_oasis` diagnostic variables to the NetCDF output | Lets you verify what wind the ocean actually received. |

**Partition/CPU‑boundary halo fixes** (received data is only valid on owned points, so it
must be halo‑filled before the physics uses stencils):

| File : line | What |
|---|---|
| `fvcom.F:1287–1339` | `ch201611 added for partition boundary issue` — a wall of `AEXCHANGE(EC,…)` on `UUWIND_oasis*`, `QPREC_cell*`, `T_AIR_cell*`, `RH_AIR_cell*`, `DLW/DSW_AIR_cell*`, `HEAT_NET_cell*`, `CLOUD_cell*`, `TDEW_cell*`, `QA_AIR_cell*`, `PA_AIR_cell*` |
| `bcond_gcn.F:654–659`, `bcond_gcy.F:542` | `ch201611 added for cpu boundary issue` — `AEXCHANGE` on `UUWIND,VVWIND,WUSURF,WVSURF` after the surface‑stress BC |
| `fvcom.F:1189` | `!!-- ch201611` — first‑two‑steps `cell0 = cell2` initialisation so the interpolation isn't garbage before the second coupling slice arrives |

**SST/AICE averaging + ice control** (`fvcom.F:986–1057`, comment `ch201611 added aice/sst
control` at 1030): SST and ice fraction are accumulated every step and divided by the period
so WRF gets the **period‑mean**, not an instantaneous snapshot. If `OASIS_ICE_CTRL_AICE`
(resp. `_SST`) is on, cells whose mean `aice ≥ OASIS_CTRL_VAL_AICE` (resp. mean
`sst ≤ OASIS_CTRL_VAL_SST`) are sent to WRF at `OASIS_freezing_temp` instead of their real
SST — a crude way to tell the atmosphere "this is ice" without a full ice‑surface scheme.

**PX‑branch merges tagged `(ch201611)`:**

| File : line | What | Role |
|---|---|---|
| `internal_step.F:1380–1435` | `px201501 … (ch201611)` deep‑water 4 °C cap, now inside `IF(DEEP_WATER_CONTROL)` (§3.4) and gated to Nov/Dec, per Great Lake, below a per‑lake depth | Keeps the hypolimnion near 4 °C in long GL runs. |
| `mod_ice.F:833–932` | `px201504 added control (ch201611)` — `IF (PRECIPITATION_ON)` around the ice precip term | Avoids using an unallocated precip array when precip forcing is off. |
| `mod_ice.F:1031–1037` | `px201504 skip ice dynamics (ch201611)` — comment out `zap_small_areas`/`rebin`/`aggregate`/`albedos`/`return` for the non‑1‑D path | Lets ice **advection/dynamics** run in the coupled GL config (the "iceDynamic" in the tree name). |
| `vertvl_edge.F:241–256` | `px201504 (ch201611)` | Vertical‑velocity edge treatment from the PX branch. |

### 3.4 `ch2017` — `NML_external_control`, runtime `ICE_MODEL`, misc init

| File : line | What | Role |
|---|---|---|
| `mod_main.F:921` | `!--ch2018 add namelist for deep-water-control` (declared here; first member `DEEP_WATER_CONTROL` is `ch2017`) — `NAMELIST /NML_external_control/` | New namelist group for Great‑Lakes‑specific controls, kept out of the standard FVCOM namelists. |
| `mod_input.F:523, 898`, `namelist.F` | `ch2017 added external control namelist` — default‑init, `READ(NML=NML_external_control)`, and `write(NML=…)` plumbing | Wires the new group into `NAME_LIST`/`NAME_LIST_INITIALIZE`/`NAME_LIST_PRINT`. |
| `mod_ncdio.F` (×4 ≈ 665, 897, 1313, 1451), `mod_report.F` (×3), `mod_startup.F` (×2), `mod_rrk.F`, `mod_enkf.F`, `mod_assim.F` (×4), `extel_edge.F`, `internal_step.F:761`, `external_step.F:105` | `IF(ICE_MODEL) THEN !--ch2017 added` | Converts the ice state (restart I/O, screen reports, startup, assimilation) from **compile‑time** `#if defined(ICE)` to the **runtime** `ICE_MODEL` logical, so one binary can run coupled with dynamic ice on or off. |
| `mod_ncdio.F:1455` | `ELSE !-ch2019 added` | The non‑ice branch of the same guard. |
| `mod_petsc.F:39,235` | `il_commlocal` (see §3.2) | Tagged `ch2017` in this file. |
| `mod_nesting.F:1547–1555` | `D%U_BLK=0.0 !ch2017 give an initial value to avoid future errors` | Zero‑init nesting block velocities (uninitialised‑memory crash in coupled runs). |
| `mod_main_wave.F:7` | `USE LIMS !ch2017` | Give the wave module the coupled dimensions. |
| `fvcom.F:802` | `sst_sum = sst_on_node !ch20170102` | Seed the period accumulators at cold start. |
| `OASIS3MCT.F_ori:120,168` | `ch2017wrf`, `ch2017` | Communicator handling in the backup coupler module. |

### 3.5 `ch201706` / `ch201706oasis` — wave → ocean radiation stress

`adv_uv_edge_gcn.F:996–1017`, `adv_uv_edge_gcy.F:800–821`, `extuv_edge.F:143–164`:

When `OASIS_wav2ocn_RadStr` is **off**, momentum flux gets the built‑in
`WAVESTRX/Y_3D` (2D) from FVCOM's own wave‑current interaction, as before. When it is **on**,
that term is suppressed and instead `field_fvcom_recv_fx/fy(i,1)*RAMP` (the radiation‑stress
gradient received from SWAN via OASIS) is added to `XFLUX/YFLUX` (internal, 3D) and
`RESX/RESY` (external, 2D). `RAMP` is FVCOM's spin‑up ramp.

### 3.6 `ch201708` / `ch201708oasis` — wave → ocean wave parameters

| File : line | What |
|---|---|
| `internal_step.F:182–243, 566, 624` | `ch201708oasis — apply coupled wave data`: under `OASIS_wav2ocn_Param`, overwrite `HSC1`, `DIRDEG1`, `TPEAK`, `WLEN`, `Pwave_bot`, `Ub_swan`, `Dwave` on owned nodes with `field_fvcom_recv_hs/dir/rtp/wlen/tmbot/ubot`, then `ACOLLECT`+`MPI_BCAST` to the global `*_GL` arrays and scatter back / `AEXCHANGE` — i.e. the coupler replaces the offline `WAVE_OFFLINE` file read. |
| `mod_main.F:1435, 2062–2071` | `ch201708oasis` — declare/allocate `hs_GL, dir_GL, rtp_GL, wlen_GL, tmbot_GL, ubot_GL` |
| `vdif_q.F:199–207` | `ch201708 added for wave offline option` — under `WAVE_OFFLINE`, take the surface wind for the TKE surface BC from `UUWIND/VVWIND` (which, when coupled, are the OASIS winds) instead of `UWWIND/VWWIND`. |

### 3.7 `ch2018` — coupler namelist reader, HRRR forcing, FVCOM 4.1 merge fixes

**Coupler‑period discovery.**

| File : line | What | Role |
|---|---|---|
| `mod_main.F:72` (`LIMS`) | `ch2018 time-variable for coupler` — `itap_sec_save`, and the six pairs `cpl_dt_<path>` / `cpl_lag_<path>` (`ocn2atm`, `ocn2wav`, `wav2ocn`, `wav2atm`, `atm2ocn`, `atm2wav`) | Per‑direction coupling period & lag, read once at init. |
| `OASIS3MCT.F:898–941` (also `.F_ori:839`) | `ch2018 — get coupler namelist info` — `subroutine get_namcpl_index(index,model_src,model_dst,cpl_dt,cpl_lag,fini)` scans OASIS's parsed `namcouple` tables (`namsrcfld`, `namdstfld`, `namflddti`, `namfldlag`, `namrstfil`) for the entry whose source starts with `model_src` and destination with `model_dst` | This is how FVCOM learns each direction's `cpl_dt`/`cpl_lag` and restart‑file name **from `namcouple` itself** rather than a duplicated FVCOM namelist. Returns `-999` / `'none'` if not found (silently disables that path — see the `namcouple` header warning). |
| `fvcom.F:535` | `ch2018` — `integer :: tmpwid,iwd,count_sst,count_aice`, `cpl_outfmt` | Loop counters / print format for the coupler block. |
| `fvcom.F:1170` | `IF (OASIS_atm2ocn) THEN !-ch2018 do only if active atm2ocn` | Guard the whole receive‑apply block. |
| `mod_par.F:252`, `mod_par_special.F:250` | `ch2018 added flag_coupler` | Coupler flag in the decomposition. |

**HRRR forcing reader** (`mod_force.F`, ~50 lines tagged `!-ch2018_hrrr`, e.g. 5680–11016):
each `FIND_DIM`/`FIND_VAR` call gets an `IF(.not. FOUND)` fallback to the HRRR/GRIB variable
name (`x`/`y`; `TMP_2maboveground`, `RH_2maboveground`, `SPFH_2maboveground`,
`PRES_surface`, `DLWRF_surface`, `DSWRF_surface`, `UGRD/VGRD_10maboveground`,
`PRATE_surface`, `TCDC_entireatmosphere`, `DPT_2maboveground` (this one `ch2019_hrrr`)…). Not
coupling per se, but it lets the **same forcing‑file code path** serve HRRR downscaling runs
and coupled runs.

**`HEATING_CALCULATED_GL`** (`mod_heatflux_gl.F:55–59`): namelist reals
`HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX`, `HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX`,
`HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX` — the GL‑tuned partition of net heat flux into
penetrating short/long wave, used when the coupler supplies radiation.
`mod_force.F:6322` `ch2018 — added for fixing seg-fault issue`: `ALLOCATE(USRCOARE(0:MT))`
for the GL COARE path.

**FVCOM 4.1 merge fixes** (behaviour parity with the group's FVCOM 3.x GL code):
`adv_uv_edge_gcn.F:67`, `adv_uv_edge_gcy.F:61`, `extuv_edge.F:55`, `mod_nctools.F:34`,
`internal_step.F:100`, `mod_ncdio.F:259` ("turn off print"), `mod_nesting.F:2255`
(single/double precision), `cntrl_prmtrs.F:173–209` (GL option + short/long‑wave
length‑scale sanity check). `chenfu2018` appears once in `mod_input.F:921` closing the
external‑control block comment.

### 3.8 `ch201808` / `ch201809` — shutdown‑hang guard

* `fvcom.F:1114, 1138` — `if (IINT<IEND) then !-ch201808` around the `fromcoupler` (RECV)
  calls for both the atm and wave paths: **don't post an OASIS receive on the final
  iteration**, because the partner model has already stopped sending → prevents a hang at
  the end of a coupled run.
* `mod_main_wave.F:36` — `!-ch201809` minor.

### 3.9 `ch2019` / `ch201902` / `ch201906` / `ch201912` — NaN guard, coupler‑namelist move, flux components

**`nan_check` (`ch201902`).** New module `nan_check.F` (`SUBROUTINE nan_check(vname)`): after
each prognostic update `internal_step.F` calls `nan_check('check'|'EL'|'UV'|'TEMP'|'SALINITY')`
(lines 133, 726, 1132, 1463, 1562); it range‑checks the field, sets `CHECK_index`, collects it
across ranks and aborts the coupled job cleanly instead of propagating NaNs into OASIS.
`namelist.F:75–263` — `!-ch201902` guards.

**Coupler namelist moved into the coupler module (`ch201902`).** `OASIS3MCT.F:1038` —
`OASIS_COUPLER_NAMELIST_INITIALIZE` / `_PRINT` / `_READ` and the `NML_OASIS_COUPLER`
definition now live in `module mod_var_cpl`, next to the code that consumes them.
`OASIS3MCT.F_ori:128` — `if (fini_fvcom(1:4)=='none') then !-ch201902` handles "no ocn→wav
entry in namcouple".

**Second wind‑speed override (`ch201902`).** `mod_force.F:13098, 13222, 16209, 16276` —
`ch201902 added for oasis wind speed`: the `WIND_TYPE == SPEED` branches of
`WINDS_ARE_WRFGRID` / `WINDS_ARE_FVCOMGRID` also get the OASIS override
(`WSTRX = UUWIND_oasis`, `WSTRXY = NEXT_WGHT_oasis*|U2| + PREV_WGHT_oasis*|U0|`), plus the
`CLOUD`/`TDEW` solar‑heating case at 13098.

**Heat‑flux component decomposition (`ch2019`, from `px201404/px2014`).**
`mod_main.F:1305, 2019` allocate `DLW, ULW, SHF, LHF`. `mod_force.F:12393–12988` — the
`SOLAR(...)` call is extended to return `tmp_nlw, tmp_shf, tmp_lhf` and
`HEAT_DLW/HEAT_ULW/HEAT_SHF/HEAT_LHF` are filled (`ch201906` heat‑flux component at 233; a
final tweak `ch202106` at 12993). `mod_ncdio.F:3504–3638` writes them (`DLW/ULW/SHF/LHF
px2014(ch2019)`). `mod_solar.F:110` — `!----ch2019`. These are diagnostics so a coupled run's
surface energy budget can be closed and compared against WRF.

**River T/S handling reverted to FVCOM 3.1.6 (`ch201906`).**
`adjust_ts.F:50–165` and `adv_t.F:597–976` — the FVCOM‑4.1 behaviour of clamping river‑node
temperature (`T1 = MAX(T1,TDIS)`) is turned **off**; salinity uses
`S1 = MAX(MIN(S1,SDIS),0)` (`ch2019`); `adv_t` takes `STPOINT = T1(JJ,K)` (upwind cell
value) "to match fvcom316". Keeps GL river inflow physics consistent across FVCOM versions.

**OBC outflow clamp disabled (`ch2019off`).** `bcond_gcn.F:351–544` (×6) — the
`UNTMP = MAX(UNTMP,0)` "no inflow at the open boundary" limiter is commented out (`!-ch2019off`),
allowing two‑way transport at the GL open boundaries. `bcond_gcn.F:686–742`,
`bcond_gcy.F:574–629` — `px201404 (ch2019 added)` OBC radiation terms.

**`USE MOD_ICE` made conditional (`ch2019 modified`).** `internal_step.F:41`,
`external_step.F:48` — `USE MOD_ICE`/`MOD_ICE2D` wrapped in `#if defined (ICE)`.

**Monthly restart (`ch201912`).** `mod_main.F:1121` `logical :: CHECK_firstday_month`;
`internal_step.F:135` sets it when the wall clock hits `01T00:00:00` (or `12-31`);
`mod_ncdio.F:286–288` forces a restart‑file write on that step. Long coupled GL runs get a
restart at the start of every month.

`mod_solar.F:110`, `swanser.F:1010` (`ch2020 turnoff warning`) round out the family.

### 3.10 `ch2020` / `ch202001` / `ch202012` — coupled restart clock

`fvcom.F:954–955`:

```fortran
! itap_sec = seconds(IMDTI)*(IINT-1)                        ! original
  itap_sec       = seconds(IMDTI)*(IINT-ISTART)             ! ch202001 solve coupled restart issue
  itap_sec_float = real(seconds(IMDTI)*(IINT-ISTART))       ! ch202012 solve the fractional dt issue (dt<1)
```

On a **hot start** the OASIS time index must restart at 0, so it is referenced to `ISTART`
(this restart's first iteration) rather than iteration 1 of the original run — otherwise
`oasis_get`/`oasis_put` look for the wrong coupling slice and the run desynchronises or
hangs. `itap_sec_float` is a real‑valued companion so that model time steps below 1 s
("frication"/fractional dt) don't get truncated to the same integer second.

### 3.11 `ch2021` / `ch202103` / `ch202104` / `ch202106` — monthly water‑level forcing

| File : line | What | Role |
|---|---|---|
| `mod_main.F:924–926` | `MONTHLY_WATER_LEVEL_RATE(12)` (`ch202103`), `MIN_DEPTH_FOR_WL_ADJUST` (`ch202104`), `DAYS_IN_MONTH(12)` (`ch202103`) added to `NML_external_control` | Prescribed lake‑level seasonal cycle (m/month) for GL runs. |
| `mod_input.F:526–528` | defaults (`0.0`, `1.0`) | — |
| `mod_main.F:1122` | `INTEGER :: index_month !ch202103` | Current month, from the wall clock. |
| `internal_step.F:104–105, 143–151, 680` | `use Control, only : MONTHLY_WATER_LEVEL_RATE, DAYS_IN_MONTH`; derive `index_month`; print the per‑day/per‑dt rate; (commented) `EL = EL + rate/DAYS_IN_MONTH/86400*seconds(IMDTI)` in the 1‑D branch | — |
| `external_step.F:356` | `ch202103 add water wlevel change rate` — the **live** application: `do i=1,M; if (h(i)>MIN_DEPTH) EL(i)=EL(i)+MONTHLY_WATER_LEVEL_RATE(index_month)/DAYS_IN_MONTH(index_month)/86400*seconds(IMDTE); end do` | Nudges free surface each external step; `MIN_DEPTH_FOR_WL_ADJUST` is a debug‑only shallow cutoff. |
| `mod_ice.F:708` | `ch2021 modified for heating calc gl` | GL heat path in the ice module. |
| `internal_step.F:127` | `real(SP) :: tmp_t1 !-ch2021` | Scratch in the deep‑water block. |

### 3.12 `ba2026` — coupling‑robustness hardening (newest)

All 2026 tags are about making the OASIS receive path fail‑safe. See §4.

---

## 4. `ba2026` in detail

The `ba2026` edits (≈25 sites in `mod_force.F`, 3 in `mod_nctools.F`) address a class of
silent‑failure bugs in the atm→ocean receive path.

### 4.1 OASIS interpolation weights were being left at 0.0 (`mod_nctools.F`)

`mod_nctools.F:9453–9480, 9606` — the block that computes `FTM%NEXT_WGHT_oasis` /
`PREV_WGHT_oasis` was **moved from the end of `UPDATE_FILE_BRACKET` to the very top**, right
after `FTM` is associated:

> *"the file‑bracket search … has several early RETURN paths (NOW outside the dummy forcing
> file's tiny bracket window, STK_LEN exhausted, etc.) that fire before reaching that code.
> When any of those paths triggered, `FTM%NEXT_WGHT_oasis`/`PREV_WGHT_oasis` were silently
> left at their compiled‑in default of 0.0/0.0 forever, which zeroed out every OASIS‑received
> variable (T_AIR, PA_AIR, wind, ice fields, etc.) for the entire run regardless of what was
> actually received from the coupler."*

Because coupled runs use tiny "dummy" forcing files (just enough to satisfy FVCOM's file
machinery), those early returns were the **normal** path — so before this fix the interpolated
fields `<V>_cell = NEXT_WGHT_oasis*<V>_cell2 + PREV_WGHT_oasis*<V>_cell0` were identically
zero. `mod_nctools.F:9439` `denom_oasis/numer_oasis` is the `ch201611` weight formula this
block still uses.

### 4.2 Skip the dummy file read, but keep the weights ticking (`mod_force.F`)

In `UPDATE_HEAT_CALCULATED`, `UPDATE_HEAT_SOLAR` (cloud/dew case), `UPDATE_WIND`
(`WRFGRID` and `FVCOMGRID`), `UPDATE_AIRPRESSURE`, and the `ICE_FILE` readers
(≈`mod_force.F:12674, 12737, 12825, 12886, 13376, 13403, 13442, 13467, 14321, 14338, 14367,
14382, 14929, 14978, 15017, 15067`):

```fortran
#if defined (oasis_coupler)
  IF ( .NOT. (USE_OASIS_COUPLER .AND. OASIS_atm2ocn) ) THEN
#endif
      ... CALL UPDATE_VAR_BRACKET(HEAT_FILE, ...)   ! file read — SKIPPED when coupled
      ... T_AIR = FTM%NEXT_WGHT*VNP + FTM%PREV_WGHT*VPP
#if defined (oasis_coupler)
  ELSE
      CALL UPDATE_FILE_BRACKET(HEAT_FILE, HTIME, STATUS)   ! weights only
  END IF
#endif
```

When coupled, the file data is never touched (the `ch201611` OASIS override just below
fully determines the field), but `UPDATE_FILE_BRACKET` is still called directly so that
§4.1's weight computation runs every timestep.

### 4.3 Skip file time‑range/dimension validation at setup

`mod_force.F:6541, 7606, 8729, 10580` — `ba2026: skip file time-range/dimension validation
entirely when OASIS-coupled`. The dummy files' time axes don't cover the run window; that's
fine because their data is unused. `HEAT_FORCING_TYPE` / grid‑type detection is still done.

### 4.4 Unconditional receive‑range diagnostics

`mod_force.F:12764, 12913, 13442` etc. — after each OASIS apply, an ungated
`print *,"OASIS RECV RANGE: T_AIR min/max", …` (and `RH_AIR`, `DLW_AIR`, `DSW_AIR`, `PA_AIR`,
`WSTRX`, `WSTRY`). Deliberately not behind `dbg_lvl`, so a bad/extreme received value is
visible in the log immediately.

### 4.5 Related but **unmarked** 2026 edits (no `ch` tag — listed for completeness)

* `fvcom.F:775–782` and `958–966` — "REORDERED" comments: FVCOM's cold‑start and
  first‑real‑step coupler blocks now **send before receive**. Both models did
  RECV‑before‑SEND, so each blocked waiting on the other at `t=0` — a permanent deadlock.
  Sending FVCOM's already‑initialised state first (SST/AICE/U/V/EL1 from `STARTUP`) lets
  WRF's blocking SST receive complete.
* `gl_run.nml` — extensive `!--` annotations explaining which `&NML_*` switches must stay
  on for the OASIS override to fire (`WIND_ON`, `HEATING_CALCULATE_ON`, `AIR_PRESSURE`,
  `PRECIPITATION_ON`), since the override lives *inside* those ON‑gated update routines.

---

## 5. Quick index by tag

| Tag | Meaning in one line | Main files |
|---|---|---|
| `ch2015`, `ch201508` | coupling time string + global error‑check scaffold | `internal_step.F`, `mod_main.F`, `mod_ncdio.F` |
| `ch201610(_oasis)` | OASIS bootstrap: communicator, `ini_coupler`, per‑step block | `fvcom.F`, `mod_main.F`(CONTROL), `mod_utils.F`, `mod_petsc.F`, `mod_par*` |
| `ch201611(_oasis)` | receive buffering, time‑interp weights, forcing override, halo fixes, SST/AICE avg + ice control, PX merges | `mod_force.F`, `mod_main.F`, `mod_nctools.F`, `mod_ncll.F`, `fvcom.F`, `bcond_gc*.F`, `mod_ice.F`, `mod_ncdio.F` |
| `ch2017` | `NML_external_control`; runtime `ICE_MODEL` guards; init fixes | `mod_input.F`, `namelist.F`, `mod_ncdio.F`, `mod_report.F`, `mod_assim.F`, `mod_nesting.F`, many |
| `ch201706(oasis)` | wave→ocean radiation‑stress path (`OASIS_wav2ocn_RadStr`) | `adv_uv_edge_gc*.F`, `extuv_edge.F` |
| `ch201708(oasis)` | wave→ocean wave‑parameter path (`OASIS_wav2ocn_Param`) | `internal_step.F`, `mod_main.F`, `vdif_q.F` |
| `ch2018` | `get_namcpl_index` (read `cpl_dt/lag` from `namcouple`); HRRR reader; `HEATING_CALCULATED_GL`; FVCOM 4.1 merge fixes | `OASIS3MCT.F`, `mod_main.F`, `mod_force.F`, `mod_heatflux_gl.F`, `cntrl_prmtrs.F` |
| `ch201808/09` | don't RECV on the last iteration (shutdown hang) | `fvcom.F`, `mod_main_wave.F` |
| `ch2019`, `ch201902` | `nan_check` module; move `NML_OASIS_COUPLER` into coupler module; 2nd wind‑speed override; heat‑flux components; `USE MOD_ICE` conditional | `nan_check.F`, `OASIS3MCT.F`, `mod_force.F`, `mod_ncdio.F`, `internal_step.F`, `external_step.F`, `namelist.F` |
| `ch201906` | river T/S handling reverted to FVCOM 3.1.6 | `adjust_ts.F`, `adv_t.F` |
| `ch2019off` | disable OBC outflow clamp `MAX(UNTMP,0)` | `bcond_gcn.F` |
| `ch201912` | write restart on the first day of each month | `mod_main.F`, `mod_ncdio.F`, `internal_step.F` |
| `ch2020`, `ch202001`, `ch202012` | coupled hot‑start clock ref'd to `ISTART`; real `itap_sec_float` for `dt<1`; SWAN warning off | `fvcom.F`, `swanser.F` |
| `ch2021`, `ch202103/04/06` | monthly prescribed water‑level forcing; GL ice heat path | `mod_main.F`, `external_step.F`, `internal_step.F`, `mod_input.F`, `mod_ice.F`, `mod_force.F` |
| `ba2026` | OASIS receive path made fail‑safe: weights computed unconditionally; skip dummy‑file reads/validation; range diagnostics | `mod_nctools.F`, `mod_force.F` |
| `chenfu2017` | `nsecrun_save`; `AICE` recv + send‑field catalogue; 2nd OASIS partition/grid; `AICE` Registry state | `frame/module_cpl*.F`, `Registry.EM_COMMON`, `module_surface_driver.F` |
| `chenfu2018` | thread `AICE` through the surface driver | `module_surface_driver.F`, `module_first_rk_step_part1.F` |
| `chenfu2019`, `chenfu201902` | per‑nest partition arrays; trimmed var names; relative winds `ru10/rv10`; `write_restart` on coupling steps; `max_cplfld=30` | `frame/module_cpl*.F`, `module_driver_constants.F` |
| `chenfu201902` (surface) | `XICE = AICE` — WRF ice fraction driven by FVCOM | `module_surface_driver.F` |
| `chenfu202005` | store SST→SST_INPUT whenever an SST field is defined | `frame/module_cpl.F` |

---