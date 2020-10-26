#!/bin/bash

nsim=$1
for (( sim = 1; sim <= nsim; sim++ )); do
	gdal_translate -of AAIGrid ../Inputs/Geostatistics/clay_sim_$sim.tif ../Inputs/Geostatistics/clay_sim_$sim.asc
	gdal_translate -of AAIGrid ../Inputs/Geostatistics/sand_sim_$sim.tif ../Inputs/Geostatistics/sand_sim_$sim.asc
done
