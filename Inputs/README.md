# Inputs/

Example input datasets. Paths in a `SMARTSED_input` file are resolved relative
to **this** directory. Gridded data is **ESRI ASCII** (`.asc`: 6-line header
`ncols/nrows/xllcorner/yllcorner/cellsize/NODATA_value`, then the raster);
`.tif` copies and `.prj`/`.aux.xml` sidecars are GIS convenience files. All
layers for one run must share a projected CRS and cover the same region.

Most of this is the real **Caldone** catchment near Lecco, Italy (UTM 32N).

| Path | Contents |
|------|----------|
| `realOrography/` | `DEM.asc` (5 m digital elevation model, 1998 × 1829) and `Mask_bin.asc` (basin mask: 1 inside, 0/NODATA outside). The reference orography for the event runs. |
| `ideal/` | Synthetic test catchment: `DEMIdeal.asc`, `mask.asc`, `CLC_id.asc`, `clay_id.asc`, `sand_id.asc`. Used by `run/SMARTSED_input_id_*`. |
| `input_smartsed/` | An alternative, self-contained input set (Cilento / Dragone–Sambuco basins): DEM, masks, CLC, soil-texture maps, rain and temperature series, Gavrilović coefficients. See its own `readme.txt`. |
| `CorineLandCover/` | `CaldoneCLC_5.asc` — CORINE Land Cover class codes on the DEM grid. Drives the SCS-CN curve numbers and the per-class Gavrilović coefficients. |
| `SoilGrids/SoilGrids2020/` | `clay.tif`, `sand.tif`, `silt.tif` — coarse particle-size-fraction maps from [SoilGrids](https://soilgrids.org). Input to the R geostatistical pre-processor, which downscales them to the DEM resolution. |
| `Geostatistics_input_data/` | The four `.tif` files (`DEM`, `clay`, `sand`, `silt`) that the R pre-processor reads — same role, drop-in location expected by the R script. |
| `DatiMeteo_ARPA_LOMBARDIA/` | Per-station meteo from ARPA Lombardia, one sub-folder per gauge (`barzio_PDB`, `canzo`, `carenno`, `cassina`, `cortenova`, `lecco`, `rotaImagna`, `sormano`, `valmadrera`). Files: `rain_<station>_ev0.txt`, `temperature_<station>_ev0.txt`. **ARPA format** = header line + `Id  YYYY/MM/DD  HH:MM  value` (rain in mm; temperature in °C). |
| `2020/` | ARPA meteo for a 2020 event, `*_bat.txt`, one per station (used by `SMARTSED_input_long_ev_2020`). |
| `temperature/`, `rain/` | Loose meteo series. `rain/event/pioggia_id.txt` and `rain/pioggia_id_constant.txt` are the two rain-file shapes (IDW multi-station vs. single uniform series). `temperature/event/temperature.txt` is the "comune" tab-delimited format (`DD/MM/YYYY  HH:MM:SS  value`). |
| `coeff_Gav.txt`, `coeff_Gav_1.txt` | Gavrilović sediment-production coefficients, two columns per CORINE land-cover class: erosion coefficient `Z` (col 1) and a land-use factor (col 2). Path is set by `physics/Gavrilovic_txt`. |

## Adding your own catchment

1. Put `DEM.tif`, `clay.tif`, `sand.tif`, `silt.tif` (same resolution/extent) in
   `Geostatistics_input_data/` and run the R pre-processor to get downscaled
   soil-texture maps — *or* provide `clay`/`sand` maps directly and set
   `restart_soilMoisture = true` (and `-sim -1`).
2. Convert the DEM and basin mask to ESRI ASCII; point `[files]` at them.
3. Provide rain and temperature series in the `arpa` or `comune` format and list
   the stations (with coordinates) in `[files/meteo_data]`.
4. Provide a CORINE land-cover raster on the DEM grid and a matching
   `coeff_Gav.txt`.
