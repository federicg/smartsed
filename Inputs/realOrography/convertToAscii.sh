#!/bin/sh

# this file convert .tif to .txt ESRI files

# the program uses AAIGrid file format


gdal_translate -of AAIGrid DEM.tif DEM.asc
gdal_translate -of AAIGrid Mask_bin.tif Mask_bin.asc

gdal_translate -of AAIGrid silt_sim_1.tif silt_sim_1.asc
gdal_translate -of AAIGrid clay_sim_1.tif clay_sim_1.asc
gdal_translate -of AAIGrid sand_sim_1.tif sand_sim_1.asc


