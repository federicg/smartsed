#!/bin/bash

nsim=$1

if [[ $nsim -ne 0 ]]; then
    for (( sim = 1; sim <= nsim; sim++ )); do
        gdal_translate -of AAIGrid ../Outputs/$sim/clay_sim_$sim.tif ../Outputs/$sim/clay_sim_$sim.asc
        gdal_translate -of AAIGrid ../Outputs/$sim/sand_sim_$sim.tif ../Outputs/$sim/sand_sim_$sim.asc
    done
else
    gdal_translate -of AAIGrid ../Outputs/$nsim/clay_sim_$nsim.tif ../Outputs/$nsim/clay_sim_$nsim.asc
    gdal_translate -of AAIGrid ../Outputs/$nsim/sand_sim_$nsim.tif ../Outputs/$nsim/sand_sim_$nsim.asc
fi
