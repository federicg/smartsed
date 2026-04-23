#!/bin/bash


NPROCS=1
JOBID=sim_ev #${PBS_JOBID}
OUTFILE=out.$NPROCS
ERRFILE=err.$NPROCS

PBS_O_WORKDIR=$PWD

rm -r ${PBS_O_WORKDIR}/../Outputs/$JOBID
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID/DeterministicProgram
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID/run
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs


cp ${PBS_O_WORKDIR}/../Outputs/*.m ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs
cp ${PBS_O_WORKDIR}/../Outputs/*.R ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs

cp -r ${PBS_O_WORKDIR}/../Outputs/initial_cond ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs/
cp -r ${PBS_O_WORKDIR}/../Geostatistics ${PBS_O_WORKDIR}/../Outputs/$JOBID/
cp -r ${PBS_O_WORKDIR}/../Inputs ${PBS_O_WORKDIR}/../Outputs/$JOBID/
cp ${PBS_O_WORKDIR}/SMARTSED_input_ev ${PBS_O_WORKDIR}/../Outputs/$JOBID/run/SMARTSED_input 
cp ${PBS_O_WORKDIR}/../DeterministicProgram/build/main.exe ${PBS_O_WORKDIR}/../Outputs/$JOBID/DeterministicProgram/main.exe


nsim=0
res=10


cd ${PBS_O_WORKDIR}/../Outputs/$JOBID/run
mpirun -np $NPROCS ../DeterministicProgram/main.exe -f SMARTSED_input -sim $nsim -scale $res >>$OUTFILE 2>$ERRFILE


#---------------------------------------------------------------------#
