#!/bin/bash
#
# Speedup table only: run CPU and CUDA builds for STEPS_SPEED accepted steps
# at several DEM resolutions and collect the [BENCH] wall-clock line to build
# a CPU-vs-GPU speedup table.
#
# Usage: ./speedup_only.sh [res1 res2 ...]

set -u

RES_LIST=("${@:-20 15 10 7 5 3 2 1}")
RES_LIST=(${RES_LIST[@]})

STEPS_SPEED=${STEPS_SPEED:-100}

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

OUTDIR=../Outputs/cpu_gpu_validate
mkdir -p "$OUTDIR"
LOG=$OUTDIR/speedup_log.txt
: > "$LOG"

echo "res list: ${RES_LIST[*]}   speed steps: $STEPS_SPEED" | tee -a "$LOG"
echo "date: $(date -Is)" | tee -a "$LOG"
echo | tee -a "$LOG"

echo "=== Speedup table ($STEPS_SPEED accepted steps) ===" | tee -a "$LOG"
printf "%-8s %-12s %-14s %-14s %-14s %-10s\n" "res" "pixel_m" "basin_cells" "cpu_wall_s" "gpu_wall_s" "speedup" | tee -a "$LOG"

for res in "${RES_LIST[@]}"; do
    STEPS=$STEPS_SPEED BUILD_DIR=build      bash submission.sh      "$res" 0 >/dev/null 2>&1
    cp "../Outputs/sim_ev_r${res}/run/out.1" "$OUTDIR/cpuN.r${res}.out"

    STEPS=$STEPS_SPEED BUILD_DIR=build-cuda bash submission_cuda.sh "$res" 0 >/dev/null 2>&1
    cp "../Outputs/sim_ev_r${res}/run/out.1" "$OUTDIR/gpuN.r${res}.out"

    cpu_wall=$(grep '^\[BENCH\]' "$OUTDIR/cpuN.r${res}.out" | sed -n 's/.*wall_s=\([0-9.]*\).*/\1/p')
    gpu_wall=$(grep '^\[BENCH\]' "$OUTDIR/gpuN.r${res}.out" | sed -n 's/.*wall_s=\([0-9.]*\).*/\1/p')
    ncells=$(grep '^\[FP\] H' "$OUTDIR/cpuN.r${res}.out" | sed -n 's/.*n=\([0-9]*\).*/\1/p')

    speedup=$(awk -v c="$cpu_wall" -v g="$gpu_wall" 'BEGIN{ if (g>0) printf "%.3f", c/g; else print "n/a" }')

    printf "%-8s %-12s %-14s %-14s %-14s %-10s\n" "$res" "$((res * 5))" "$ncells" "$cpu_wall" "$gpu_wall" "$speedup" | tee -a "$LOG"
done

echo | tee -a "$LOG"
echo "done -> $LOG" | tee -a "$LOG"
