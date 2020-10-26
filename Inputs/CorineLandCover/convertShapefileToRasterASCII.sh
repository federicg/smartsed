#!/bin/bash


scale_factor=$1
res_final=$(($scale_factor*5))
for (( res = 5; res <= $res_final; res=res+5 )); do
	echo "$res"
	gdal_rasterize -l CaldoneCLC -a CLC12_3L_2 -tr $res $res -a_nodata 999999.0 -te 528669.6465999996 5076881.926100001 538659.6465999996 5086026.926100001 -ot Float32 -of GTiff ../Inputs/CorineLandCover/CaldoneCLC.shp ../Inputs/CorineLandCover/CaldoneCLC_$res.tif
	gdal_translate -of AAIGrid ../Inputs/CorineLandCover/CaldoneCLC_$res.tif ../Inputs/CorineLandCover/CaldoneCLC_$res.asc
done


