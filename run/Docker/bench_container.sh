#!/bin/bash
#
# Entry point of the smartsed-cuda image: runs the event benchmark on the CPU
# (1 core, Release build) and on the GPU, and reports the speed-up.
#
# Usage:  smartsed-bench [res ...]          default: 20 10 5 3 2 1
#   res   DEM coarsening factor, pixel = res * 5 m (res 1 = native 5 m grid)
# Env:
#   STEPS  accepted time steps per run (default 100)
#   STEPS_R1  accepted time steps at res=1, where the GPU loop is short (default 200)
#   CPU    set to 0 to skip the CPU baseline runs
#   WALL   if set (seconds), run res=1 without a step limit for WALL seconds
#          (includes ~15 min of setup at 5 m, so use e.g. WALL=3600)
#          and report simulated-time / wall-time (the 10-day estimate)
#
# Everything is written under /results (bind-mount a host directory there):
#   /results/summary.txt   one block per run: GPU, [BENCH] line, peak GPU memory
#   /results/work/         staged inputs + per-run logs

set -u
RES_LIST=(${@:-20 10 5 3 2 1})
STEPS=${STEPS:-100}
STEPS_R1=${STEPS_R1:-200}
SRC=${SRC:-/opt/smartsed}
OUT=${OUT:-/results}
W=$OUT/work
SUM=$OUT/summary.txt

# main.exe resolves ../Inputs, ../Outputs, ../DeterministicProgram from its run dir
mkdir -p $W/Outputs
cp -r --update=none $SRC/Inputs $SRC/DeterministicProgram $W/
cp -r --update=none $SRC/Outputs/initial_cond $W/Outputs/

{
  echo "date: $(date -Is)"
  nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version --format=csv,noheader
  echo "steps: $STEPS (res 1: $STEPS_R1)  res: ${RES_LIST[*]}  wall: ${WALL:-none}"
  echo
} | tee $SUM

run_one () {   # $1 = CPU|GPU, $2 = res, rest = extra args
  local tag=$1 res=$2; shift 2
  local exe=main.exe; [ $tag = CPU ] && exe=main_cpu.exe
  local d=$W/run_${tag}_r$res
  rm -rf $d; mkdir -p $d; cd $d
  # gauge time series are written every step by the CPU build only: keep file I/O out of the timed loop
  sed 's/^\( *save_temporal_sequence *=\) *true/\1 false/' $SRC/run/SMARTSED_input_ev > SMARTSED_input
  ../DeterministicProgram/$exe -f SMARTSED_input -sim 0 -scale $res "$@" > out 2> err &
  local pid=$! peak=0 m
  while kill -0 $pid 2>/dev/null; do
    m=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    [ "${m:-0}" -gt $peak ] && peak=$m
    sleep 0.5
  done
  wait $pid; local rc=$?
  echo "=== $tag res=$res  pixel=$((res * 5)) m  exit=$rc  peak_gpu_mem=${peak} MiB (whole GPU) ===" | tee -a $SUM
  grep -E '^\[BENCH\]|^\[FP\] H ' out | tee -a $SUM
  [ $rc -ne 0 ] && tail -n 5 err | tee -a $SUM
  echo | tee -a $SUM
}

wall () { grep -o 'wall_s=[0-9.]*' $W/run_$1_r$2/out 2>/dev/null | cut -d= -f2; }

for res in "${RES_LIST[@]}"; do
  n=$STEPS; [ $res = 1 ] && n=$STEPS_R1
  [ "${CPU:-1}" != 0 ] && run_one CPU $res -steps $n
  run_one GPU $res -steps $n
  c=$(wall CPU $res) g=$(wall GPU $res)
  [ -n "$c" ] && [ -n "$g" ] && \
    awk -v r=$res -v c=$c -v g=$g 'BEGIN{printf "=== speed-up res=%s: CPU %.2f s / GPU %.2f s = %.2fx ===\n\n", r, c, g, c/g}' | tee -a $SUM
done

if [ -n "${WALL:-}" ]; then
  # max_Days in SMARTSED_input_ev is 1, so "Simulation progress" is % of one day
  d=$W/run_r1_wall
  rm -rf $d; mkdir -p $d; cd $d
  cp $SRC/run/SMARTSED_input_ev SMARTSED_input
  timeout $WALL ../DeterministicProgram/main.exe -f SMARTSED_input -sim 0 -scale 1 2> err \
    | awk '/Simulation progress/{print systime(), $3; fflush()}' > progress.txt
  # rate over the second half of the loop only: the start-up transient has tiny time steps
  awk '{t[NR]=$1; p[NR]=$2}
       END{h=int(NR/2); if (h < 1 || t[NR] <= t[h]) {print "=== res=1 timed run: no time-loop progress, raise WALL ==="; exit}
           r=(p[NR]-p[h])*864/(t[NR]-t[h]);
           printf "=== res=1 timed run: loop %d s, simulated %.2f%% of a day, %.2fx real time (2nd half), 10 days ~ %.1f h ===\n",
           t[NR]-t[1], p[NR], r, 864000/r/3600}' progress.txt | tee -a $SUM
fi

echo "done -> $SUM"
