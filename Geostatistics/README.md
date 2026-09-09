# Geostatistics/ — soil-texture pre-processor

Produces the downscaled clay / silt / sand maps that the deterministic solver
needs, starting from coarse [SoilGrids](https://soilgrids.org) particle-size
fractions and the catchment DEM. Geostatistical downscaling is done in the
**Aitchison geometry** (the fractions are compositional data — they sum to 1),
by kriging with an external drift on the DEM, with optional conditional
simulations for uncertainty.

## Files

| File | Role |
|------|------|
| `pullSoilGridsdata.py` | Downloads SoilGrids raster tiles (clay/sand/silt, 0–5 cm mean) via the ISRIC WCS service, using the **QGIS** Python API. Needs QGIS installed and `PYTHONPATH` / `PATH` pointed at it (see the header comment). Paths inside are hard-coded — edit before use. |
| `DownscalingAitchisonSmartSed_2020.R` | The downscaling itself. Reads `../Inputs/Geostatistics_input_data/DEM.tif` as the target grid and the coarse SoilGrids maps, transforms to Aitchison coordinates, krigs/simulates onto the DEM resolution, back-transforms, and writes the fine clay/sand/silt maps. Takes the number of conditional simulations as `args[1]`. |

## Inputs expected

In `../Inputs/Geostatistics_input_data/`:

- `DEM.tif` — digital elevation model; defines the target resolution, extent, CRS.
- `clay.tif`, `sand.tif`, `silt.tif` — coarse SoilGrids fractions, same
  resolution as each other, covering the same region as the DEM.

(`../Inputs/SoilGrids/SoilGrids2020/` holds a copy of the coarse maps as
downloaded.)

## Dependencies

- **R**: `raster`, `gstat`, `compositions`, `dissever`, `fields`, `soiltexture`,
  `viridis`, `psych`, `latex2exp`
- **GDAL** binaries on the system
- **QGIS** (only for `pullSoilGridsdata.py`; the download can also be done by
  hand from <https://soilgrids.org> or <https://files.isric.org/soilgrids/data/recent/>)

## Where it fits

```
SoilGrids WCS ──pullSoilGridsdata.py──▶ Inputs/Geostatistics_input_data/*.tif
                                              │
                     DownscalingAitchisonSmartSed_2020.R
                                              │
                                              ▼
                     fine clay/sand/silt maps ──▶ SMARTSED_input
                                                  [files/initial_conditions] clay_file / sand_file
```

The solver's `run/submission*.sh` scripts reference a `../Geostatistics`
directory relative to `run/` — i.e. this folder. Run this step first (or set
`restart_soilMoisture = true` with `-sim -1` to skip it and supply the maps
directly).
