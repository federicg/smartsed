#!/bin/bash
#
# CPU reference run (build-debug, profiled with -pg).
#
# Usage:  ./submission.sh [res] [nsim]
#   res   DEM coarsening factor passed to -scale.  pixel_size = res * 5 m
#         (base DEM cellsize is 5 m, 1998 x 1829).  Default: 10  (-> 50 m grid)
#   nsim  Monte Carlo simulation index passed to -sim.  Default: 0
#
# Env:
#   STEPS      if set (>0), append "-steps $STEPS" so the run stops after that
#              many accepted time steps (benchmark mode).  Default: full run.
#   BUILD_DIR  DeterministicProgram build dir to deploy.  Default: build-debug
#
# Output goes to  ../Outputs/$JOBID  with JOBID = sim_ev_r<res>  so that runs
# at different resolutions do not overwrite each other (needed for the
# CPU-vs-GPU speedup sweep).

NPROCS=1

res=${1:-10}
nsim=${2:-0}

JOBID=sim_ev_r${res}
PBS_O_WORKDIR=$PWD

RUNDIR=${PBS_O_WORKDIR}/../Outputs/$JOBID/run
OUTFILE=$RUNDIR/out.$NPROCS
ERRFILE=$RUNDIR/err.$NPROCS

rm -rf ${PBS_O_WORKDIR}/../Outputs/$JOBID
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID/DeterministicProgram
mkdir -p $RUNDIR
mkdir -p ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs

cp ${PBS_O_WORKDIR}/../Outputs/*.m ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs 2>/dev/null
cp ${PBS_O_WORKDIR}/../Outputs/*.R ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs 2>/dev/null
cp -r ${PBS_O_WORKDIR}/../Outputs/initial_cond ${PBS_O_WORKDIR}/../Outputs/$JOBID/Outputs/
cp -r ${PBS_O_WORKDIR}/../Geostatistics ${PBS_O_WORKDIR}/../Outputs/$JOBID/
cp -r ${PBS_O_WORKDIR}/../Inputs ${PBS_O_WORKDIR}/../Outputs/$JOBID/
cp ${PBS_O_WORKDIR}/SMARTSED_input_ev $RUNDIR/SMARTSED_input

# BUILD_DIR can be overridden (e.g. build / build-debug / build-cuda)
BUILD_DIR=${BUILD_DIR:-build-debug}
cp ${PBS_O_WORKDIR}/../DeterministicProgram/${BUILD_DIR}/main.exe \
   ${PBS_O_WORKDIR}/../Outputs/$JOBID/DeterministicProgram/main.exe

# Activate conda
source ~/HDD/miniforge3/etc/profile.d/conda.sh
conda activate base

cd $RUNDIR

# Set this env var to force gmon.out location
export GMON_OUT_PREFIX=$RUNDIR/gmon.out

STEPS_ARG=""
[ -n "${STEPS:-}" ] && [ "${STEPS:-0}" -gt 0 ] 2>/dev/null && STEPS_ARG="-steps $STEPS"

echo "[submission] res=$res  nsim=$nsim  pixel_size=$((res * 5)) m  ${STEPS_ARG:-full-run}" | tee -a $OUTFILE

mpirun -np $NPROCS \
    ../DeterministicProgram/main.exe -f SMARTSED_input -sim $nsim -scale $res $STEPS_ARG \
    >>$OUTFILE 2>$ERRFILE

echo "Exit code: $?"
ls -lah $RUNDIR/gmon.out* 2>/dev/null || echo "No gmon.out found"
