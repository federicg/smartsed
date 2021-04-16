#!/bin/bash


if [[ $# -ne 2 ]]; then
	echo 'Provide two and only two input arguments!'
    exit 1
fi


nsim=$1 # number of stochastic simulations of zero only kriging performed
res=$2

# Sequential simulations
# Geostatistics:
#mkdir -p ../Inputs/Geostatistics
#OUTFILE=../Inputs/Geostatistics/outStochastic.txt
#ERRFILE=../Inputs/Geostatistics/errStochastic.txt

if [[ $nsim -ge 0 ]]; then
  Rscript ../Geostatistics/Downscaling_Simulation_SoilGrids/Downscaling/DownscalingAitchisonSmartSed_2020.R $nsim $res
fi

chmod +x ../Geostatistics/Downscaling_Simulation_SoilGrids/Downscaling/convertToRasterASCII.sh
./../Geostatistics/Downscaling_Simulation_SoilGrids/Downscaling/convertToRasterASCII.sh $nsim

chmod +x ../Inputs/CorineLandCover/convertShapefileToRasterASCII.sh
./../Inputs/CorineLandCover/convertShapefileToRasterASCII.sh $res




# Deterministic part: 
if [[ $nsim -gt 0 ]]; then
	for (( i = 1; i <= $nsim; i++ )); do
		OUTFILE=../Outputs/$i/outdeterministic.txt
		ERRFILE=../Outputs/$i/errDeterministic.txt
	    OMP_NUM_THREADS=4 ./../DeterministicProgram/main.exe -f SMARTSED_input -sim $i -scale $res #>>$OUTFILE 2>$ERRFILE
	done
else
	OUTFILE=../Outputs/$nsim/outdeterministic.txt
	ERRFILE=../Outputs/$nsim/errDeterministic.txt
    OMP_NUM_THREADS=1 ./../DeterministicProgram/main.exe -f SMARTSED_input -sim $nsim -scale $res #>>$OUTFILE 2>$ERRFILE
fi

