# SMART-SED

**SMART-SED** is a distributed, physically-based model of the **hydrological and
erosion/sediment-transport response of a mountain catchment**. Given a digital
elevation model, a basin mask, soil-texture maps and meteorological time series
(rain, temperature), it integrates the 2D shallow-water (de Saint-Venant)
equations coupled to snow accumulation/melt, infiltration, evapotranspiration, a
sub-surface gravitational water layer and Gavrilović-type sediment production,
and writes gridded maps and gauge time series of water depth, velocity, snow,
soil moisture and sediment.

The project has two parts:

| Part | Language | What it does |
|------|----------|--------------|
| **Geostatistical pre-processor** | R | Downscales coarse SoilGrids clay/silt/sand maps to the DEM resolution (kriging + optional conditional simulations) to produce the soil-texture inputs. |
| **Deterministic solver** | C++17 (+ optional CUDA), MPI | Time-marches the coupled PDE system on the basin grid. This is the code in [`DeterministicProgram/`](DeterministicProgram/). |

MPI is used only to run several **stochastic realisations** in parallel (one rank
per realisation); each single simulation is serial (or GPU-accelerated).

For the full modelling background and a parameter-by-parameter walk-through of the
input file, see **[`Userguide_SMARTSED.pdf`](Userguide_SMARTSED.pdf)**.

---

## Repository layout

```
smartsed/
├── DeterministicProgram/   C++/CUDA solver + bundled headers (Eigen IML++, GetPot)
│                           see DeterministicProgram/README.md
├── Inputs/                 example input datasets (DEM, masks, meteo, soil, land cover)
│                           see Inputs/README.md
├── Outputs/                simulation output goes here at run time (git-ignored);
│                           Outputs/initial_cond/ holds checked-in restart fields
│                           see Outputs/README.md
├── run/                    run configuration + launch scripts + Docker environment
│                           see run/README.md
├── report/                 LaTeX report on the CUDA port + benchmark
│                           see report/README.md
├── Userguide_SMARTSED.pdf  end-user manual (modelling + input reference)
└── README.md               this file
```

`Zeus/` and `contratto/` (if present in your working copy) are local, private
material and are git-ignored — do not commit them.

---

## Dependencies

### Deterministic solver (C++)
- CMake ≥ 3.18, a C++17 compiler
- **MPI** (OpenMPI / MPICH)
- **Eigen 3** (≥ 3.3, `find_package(Eigen3 ... NO_MODULE)`)
- *optional GPU:* CUDA Toolkit ≥ 12 (cuSPARSE, cuBLAS), an NVIDIA GPU
  (default target: `sm_86` / Ampere — edit `CMAKE_CUDA_ARCHITECTURES` in
  `DeterministicProgram/CMakeLists.txt` for other cards), and **NetCDF-C**
  (`nc-config` on `PATH`; the GPU path writes NetCDF instead of ESRI ASCII).

Eigen's IML++ Krylov solvers and the GetPot parser are vendored under
`DeterministicProgram/include/` — nothing to install for those.

### Geostatistical pre-processor (R)
R with: `raster`, `gstat`, `compositions`, `dissever`, `fields`, `soiltexture`,
`viridis`, `psych`, `latex2exp`, plus the **GDAL** binaries on the system.

### Docker (recommended on Windows / for a reproducible environment)
A `Dockerfile` that installs all of the above (minus CUDA) lives in
[`run/Docker/`](run/Docker/).

---

## Build (deterministic solver)

All commands are run from the repository root and target
`DeterministicProgram/` as the CMake source dir.

```sh
# Release (optimised, CPU) — the default
cmake -S DeterministicProgram -B DeterministicProgram/build
cmake --build DeterministicProgram/build -j

# Debug (CPU, -O0 -g, profiling with -pg)
cmake -S DeterministicProgram -B DeterministicProgram/build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build DeterministicProgram/build-debug -j

# Release + CUDA
cmake -S DeterministicProgram -B DeterministicProgram/build-cuda -DENABLE_CUDA=ON
cmake --build DeterministicProgram/build-cuda -j

# Debug + CUDA
cmake -S DeterministicProgram -B DeterministicProgram/build-cuda-debug \
      -DENABLE_CUDA=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build DeterministicProgram/build-cuda-debug -j
```

Each build produces `main.exe` in its build directory. The same source tree
compiles with and without `-DENABLE_CUDA` — the CUDA and CPU code paths are
selected at compile time with `#ifdef ENABLE_CUDA`, so the CPU build is an exact
reference for validating the GPU build.

---

## Run

The solver is always invoked through a run directory that contains a
`SMARTSED_input` configuration file. The helper scripts in `run/` set that up
for you:

```sh
cd run
./submission.sh        [res] [nsim]    # CPU build   (build-debug by default)
./submission_cuda.sh   [res] [nsim]    # CUDA build  (build-cuda   by default)
```

- `res`  – DEM coarsening factor; grid spacing = `res × 5 m` (base DEM is 5 m).
- `nsim` – stochastic-realisation index passed to `-sim`.

Output is written to `Outputs/sim_ev_r<res>/`. See [`run/README.md`](run/README.md)
for the input-file reference, the direct `mpirun` invocation, the CLI flags
(`-f`, `-sim`, `-scale`, `-steps`) and the CPU-vs-GPU benchmark harness.

---

## Authors

Federico Gatti — `federico.gatti@math.ethz.ch` (deterministic solver, CUDA port).
SMART-SED was developed at MOX, Politecnico di Milano.

## Code formatting

C++ is formatted with `clang-format` (LLVM-ish style). From vim: `:!clang-format -i %`.
