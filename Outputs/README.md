# Outputs/

Simulation results are written here at run time. The `run/submission*.sh`
scripts create one sub-directory per run, `Outputs/sim_ev_r<res>/` (and
`Outputs/bench_sweep/` for the benchmark harness). **These are git-ignored** —
only `initial_cond/` and this file are tracked.

## `initial_cond/`

Checked-in restart fields for the reference event (`SMARTSED_input_ev` points
its `[files/initial_conditions]` here):

| File | Field |
|------|-------|
| `H_ev.asc`, `h_ev_5.asc` | water surface depth (m), coarse / 5 m |
| `u_ev.asc`, `v_ev.asc`, `u_ev_5.asc`, `v_ev_5.asc` | depth-averaged velocity components (m/s) |
| `hG_ev.asc` | gravitational (sub-surface) water-layer depth (m) |
| `hsd_ev.asc` | sediment accumulation (m, normalised to cell size) |
| `clay_0.asc`, `sand_0.asc` | soil particle-size fractions (downscaled) |
| `elab.R` | R script that post-processes / resamples these fields and makes the water-height and velocity plots (`raster`, `ggquiver`). |

## Output of a run

**CPU build** — ESRI ASCII maps every `debug/frequency_save` hours, with an
index suffix:

| Prefix | Field (units) |
|--------|---------------|
| `H_<i>` | water surface depth (m) |
| `u_<i>`, `v_<i>` | velocity components (m/s), on the staggered grid |
| `hsd_<i>` | sediment accumulation (m, per cell) |
| `hG_<i>` | gravitational water-layer depth (m) |
| `hsn_<i>` | snow depth (m) |
| `w_cum_<i>` | cumulative sediment-source map (m) |

Plus, if `discretization/save_temporal_sequence = true`, two-column
`<name>.txt` gauge time series (**value, time-in-seconds**) at each
`X_gauges_i`/`Y_gauges_i`: `waterSurfaceHeight` (m), `waterSurfaceMassFlux`
(m²/s), `SolidFlux` (m²/s), `timesteps` (s).

**CUDA build** — the same fields written to **NetCDF** (`saveToNetCDF`) instead
of per-field ASCII.

Post-process the ASCII maps in QGIS / MATLAB; the time series in Python / R / MATLAB.
