# DeterministicProgram — the SMART-SED solver

C++17 (+ optional CUDA), MPI. Time-marches the coupled shallow-water /
snow / infiltration / gravitational-layer / sediment system on the basin grid
and writes gridded maps and gauge time series.

Build instructions are in the [top-level README](../README.md#build-deterministic-solver);
run instructions are in [`run/README.md`](../run/README.md).

---

## Source files

| File | Role |
|------|------|
| `main_final_H.cpp` | `main()`. Reads the input file (GetPot), builds the mesh / adjacency / staggered-grid index sets, sets up initial and boundary conditions, then runs the time-stepping loop: fluxes → friction → meteo (rain/temperature/ET) → infiltration → snow & gravitational layer → velocity interpolation → sediment sub-loop → assemble & solve the DSV linear system → update `H`, `η`, `u`, `v` → adaptive Δt. Writes output at `frequency_save` intervals. |
| `code_init.h` / `code_init.cpp` | Grid & topology setup used before the time loop: `Vector2D`, `computeAdjacencies`, `computePourCell`, pour-point / drainage logic, sink filling, ESRI-ASCII raster readers, `saveSolution` (writes an output map), `compute_d_perc` (soil particle-size fractions). Pure CPU, no CUDA. |
| `utils_H.h` / `utils_H.cpp` | The physics kernels as templated CPU functions **and** the host-side declarations of their GPU wrappers. Key pieces: `class Rain` (IDW / constant precipitation + SCS-CN infiltration), `class Temperature` (lapse-rate temperature, Hargreaves ET), `computeResiduals` / `computeResidualsTruncated` (gravitational-layer and sediment fluxes), `bilinearInterpolation` (cell→staggered velocity), `buildMatrix` (DSV sparse system assembly), `updateVel`, `compute_dt_adaptive`, `maxdt`, `compute_dt_sediment`. Also defines Matrix-Market I/O helpers for debugging. |
| `cuda_utils_loop_H.cu` / `.cuh` | The CUDA port: ~38 `__global__` kernels behind ~17 host wrappers (`*_wrapper`), device reductions (`deviceMax/Min/Sum`), the cuSPARSE/cuBLAS preconditioned-CG solver with a red-black IC(0) preconditioner, and NetCDF output (`saveToNetCDF`). Only compiled/linked when `-DENABLE_CUDA=ON`. |
| `timing.h` | Minimal scoped timers, active when compiled with `-DENABLE_TIMING` (always on — see `CMakeLists.txt`). |
| `CMakeLists.txt` | Single project, `option(ENABLE_CUDA)`. When CUDA is on, the `.cpp` files are recompiled with `nvcc`, and NetCDF is located via `nc-config`. |
| `include/` | Vendored third-party headers — nothing to install. |

### `include/`

| Path | What |
|------|------|
| `include/GetPot.hpp` | Header-only command-line / config-file parser. Parses `SMARTSED_input`. |
| `include/eigen_forma/*.hpp` | IML++ Krylov solver templates (`cg`, `bicgstab`, `gmres`, `minres`, …) adapted to Eigen. Only `cg.hpp` is used by the CPU build. |
| `include/eigen_forma/Utilities/` | COO matrix extractor, Matrix-Market readers, `RotatingVector`, graph/reordering (RCM) utilities, with their own `README.md` and `test/`. Used for debugging the assembled linear system. |

---

## CPU vs. GPU

The same translation units compile both ways:

- **`-DENABLE_CUDA=OFF`** (default) — everything runs on the host; the assembled
  DSV system is solved with Eigen + IML++ CG (or a direct `SimplicialLDLT` if
  `linear_solver/direct_method = true`). Output is ESRI ASCII (`.asc`).
- **`-DENABLE_CUDA=ON`** — the entire inner time loop runs on the device
  (interface fluxes, friction, meteo, `computeResiduals`, snow/grav update,
  velocity interpolation, sediment sub-loop, sparse assembly, PCG solve, `H`/`u`/`v`
  update, adaptive Δt). Output is NetCDF.

The CPU build is the numerical reference for the GPU build. A field-fingerprint
harness (the `[FP]` lines printed with `-steps N`) compares the two backends
step by step; with `--use_fast_math` they agree to ~1e-7 per step rather than
bit-for-bit. See [`../report/`](../report/) for the full port status, validation
and benchmark.

---

## Command-line flags

`main.exe` is launched via `mpirun` (see `run/`). Flags (parsed by GetPot):

| Flag | Meaning | Default |
|------|---------|---------|
| `-f <file>` / `--file <file>` | path to the `SMARTSED_input` config file | `SMARTSED_input` |
| `-sim <n>` | stochastic-realisation index (also selects MPI work split) | `2` |
| `-scale <k>` | DEM coarsening factor; grid spacing = `k × 5 m` | `2` |
| `-steps <N>` | **benchmark mode:** stop after `N` accepted time steps and print `[BENCH]` wall-clock + `[FP]` field fingerprints (overrides `debug/max_steps`) | `0` = full run |

Additional debug switches live in the `[debug]` block of the input file
(`stop_after_first_step`, `spit_out_matrix`, `spit_out_solutions_each_time_step`, …).
