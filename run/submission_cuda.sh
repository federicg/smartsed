#!/bin/bash

NPROCS=1
JOBID=sim_ev
PBS_O_WORKDIR=$PWD

RUNDIR=${PBS_O_WORKDIR}/../Outputs/$JOBID/run
OUTFILE=$RUNDIR/out.$NPROCS
ERRFILE=$RUNDIR/err.$NPROCS

rm -r ${PBS_O_WORKDIR}/../Outputs/$JOBID
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID/DeterministicProgram
mkdir -p $RUNDIR
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs

cp ${PBS_O_WORKDIR}/../Outputs/*.m ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs
cp ${PBS_O_WORKDIR}/../Outputs/*.R ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs
cp -r ${PBS_O_WORKDIR}/../Outputs/initial_cond ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs/
cp -r ${PBS_O_WORKDIR}/../Geostatistics ${PBS_O_WORKDIR}/../Outputs/$JOBID/
cp -r ${PBS_O_WORKDIR}/../Inputs ${PBS_O_WORKDIR}/../Outputs/$JOBID/
cp ${PBS_O_WORKDIR}/SMARTSED_input_ev $RUNDIR/SMARTSED_input
cp ${PBS_O_WORKDIR}/../DeterministicProgram/build-cuda/main.exe \
   ${PBS_O_WORKDIR}/../Outputs/$JOBID/DeterministicProgram/main.exe

nsim=0
res=10

# Activate conda
source ~/HDD/miniforge3/etc/profile.d/conda.sh
conda activate base

cd $RUNDIR

mpirun -np $NPROCS \
    ../DeterministicProgram/main.exe -f SMARTSED_input -sim $nsim -scale $res \
    >>$OUTFILE 2>$ERRFILE
