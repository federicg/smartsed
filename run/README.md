# run/ — configuration and launch

Everything needed to launch a SMART-SED simulation: the input files that
configure a run, the scripts that stage a run directory and call `mpirun`, the
benchmark harness, and a Docker environment.

---

## Contents

| File | Purpose |
|------|---------|
| `SMARTSED_input_ev` | Reference **event** run (~1 day, real Caldone basin, ARPA meteo). Used by `submission*.sh` and the benchmark. |
| `SMARTSED_input_id_1`, `SMARTSED_input_id_2` | **Ideal** test cases (synthetic DEM in `Inputs/ideal/`), sediment transport off, long horizon. |
| `SMARTSED_input_long_ev_2020` | Long 2020 event, ARPA Lombardia meteo (`Inputs/2020/`). |
| `SMARTSED_input_year` | Full-year run (`max_Days = 365`), sediment on. |
| `SMARTSED_review` | Short run (`max_Days = 0.125`) used for review/debugging. |
| `submission.sh` | Stage + launch the **CPU** build (`build-debug` by default, profiled with `-pg`). |
| `submission_cuda.sh` | Stage + launch the **CUDA** build (`build-cuda` by default). |
| `benchmark_sweep.sh` | Run both builds over a list of resolutions and collect the `[BENCH]`/`[FP]` lines into `../Outputs/bench_sweep/summary.txt`. |
| `Docker/` | Ubuntu 20.04 image with all C++ and R dependencies (no CUDA). `buildDocker*` builds it, `runUnix.sh` / `runWINDOWS.bat` mount the repo at `/smartsed` and open a shell. |

---

## Running

### Via the helper scripts (recommended)

```sh
cd run
./submission.sh       [res] [nsim]     # CPU
./submission_cuda.sh  [res] [nsim]     # GPU
```

- **`res`** — DEM coarsening factor. Grid spacing = `res × 5 m` (base DEM cell
  size is 5 m, 1998 × 1829). Default `10` → 50 m grid.
- **`nsim`** — stochastic-realisation index passed to `-sim`. Default `0`.

Environment variables honoured by both scripts:

| Var | Effect |
|-----|--------|
| `STEPS` | if `> 0`, append `-steps $STEPS` → stop after that many accepted time steps and print benchmark lines. |
| `BUILD_DIR` | which `DeterministicProgram/<build>/` tree to deploy (`build`, `build-debug`, `build-cuda`, …). |

Each script wipes and recreates `../Outputs/sim_ev_r<res>/`, copies `Inputs/`,
`Outputs/initial_cond/` and the chosen `SMARTSED_input_*` into it, copies
`main.exe` from the build dir, then runs. Results land in
`../Outputs/sim_ev_r<res>/run/` (`out.1`, `err.1`).

> The scripts `source` a conda profile (`~/HDD/miniforge3/...`) and reference a
> `../Geostatistics` directory — both are environment-specific; edit the paths
> for your machine, or run `main.exe` directly (below).

### Directly

```sh
mkdir myrun && cd myrun
cp ../run/SMARTSED_input_ev SMARTSED_input
mpirun -np 1 /path/to/DeterministicProgram/build/main.exe \
       -f SMARTSED_input -sim 0 -scale 10
```

`-np N` runs `N` stochastic realisations in parallel (one per rank); each
realisation is itself serial (or GPU). Paths inside the input file are resolved
relative to `../Inputs/`, so the run directory must sit one level below a copy
of `Inputs/` (this is what the helper scripts arrange).

### Benchmark sweep

```sh
cd run
STEPS=100 ./benchmark_sweep.sh 20 10 5 3      # res list; pixel = res × 5 m
# -> ../Outputs/bench_sweep/summary.txt
```

---

## The `SMARTSED_input` file

GetPot / INI-style: `[section]` … `[../]` nesting, `#` comments. Paths are
relative to `Inputs/`. Full reference with units and semantics is in
[`../Userguide_SMARTSED.pdf`](../Userguide_SMARTSED.pdf); summary:

| Section | Key parameters |
|---------|----------------|
| `[files]` | `orography_file`, `mask_file` (both ESRI ASCII). |
| `[files/meteo_data]` | `temperature_file`, `height_thermometer`, `format_temp` (`arpa`\|`comune`); `precipitation`, `constant_precipitation` (false ⇒ IDW), `number_stations` and per-station `rain_file_i` / `X_i` / `Y_i` / `time_spacing_rain_i`. |
| `[files/initial_conditions]` | `restart_{H,vel,snow,sediment,gravitational,soilMoisture}` flags + the corresponding `*_file` maps (finer-resolution maps are auto-resampled). |
| `[files/infiltration]` | `infiltration_model` (`SCS-CN`\|`None`), `isInitialLoss`, `perc_initialLoss`, `roughness_scale_factor{1,2,3}` (by slope class), `corineCode_file`. |
| `[files/evapotranspiration]` | `ET_model` (`Hargreaves`\|`None`), `latitude_deg`. |
| `[physics]` | `friction_model` (`Manning`\|`Rickenmann`\|`None`), `n_manning`, `Gavrilovic_txt`, `is_sediment_transport`. |
| `[discretization]` | `FillSinks`, `steps_per_hour`, `max_Days`, `starting_day`, `H_min`, `T_thr` (snow/rain threshold °C), `number_gauges` + `X_gauges_i`/`Y_gauges_i`, `delta_gauges`, `save_temporal_sequence`, `isNonReflectingBC`, `slope_thr`, `static_subbasin_approx`. |
| `[linear_solver]` | `direct_method` (true ⇒ Cholesky, false ⇒ preconditioned CG), `use_preconditioner`. |
| `[debug]` | `frequency_save` (hours between output maps), `spit_out_solutions_each_time_step`, `stop_after_first_step`, `max_steps`, `spit_out_matrix`, `matrix_name`, `vector_name`. |
