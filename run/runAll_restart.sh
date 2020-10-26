#!/bin/bash


if [[ $# -ne 2 ]]; then
	echo 'Provide two and only two input argument!'
    exit 1
fi


nsim=$1 # number of stochastic simulations
res=$2






#docker run --rm -it -v `pwd`:/smartsed smartsed_env make -f Makefile_windows

# Deterministic part: 
for (( i = 1; i <= $nsim; i++ )); do
	mkdir -p ../Outputs/$i
	OUTFILE=../Outputs/$i/outdeterministicRestart.txt
	ERRFILE=../Outputs/$i/errDeterministicRestart.txt
    OMP_NUM_THREADS=4 ./../DeterministicProgram/main.exe -f SMARTSED_input_restart -sim $i -scale $res #>>$OUTFILE 2>$ERRFILE
done


