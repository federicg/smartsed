#!/bin/bash
#
# 1) Correctness check: run CPU and CUDA builds for exactly one accepted time
#    step (STEPS=1) at several DEM resolutions and diff the [FP] field
#    fingerprint lines main_final_H.cpp prints at the end of the (truncated)
#    time loop.
# 2) Speedup table: run CPU and CUDA builds for STEPS_SPEED accepted steps at
#    the same resolutions and collect the [BENCH] wall-clock line to build a
#    CPU-vs-GPU speedup table.
#
# Usage: ./cpu_gpu_validate_and_speedup.sh [res1 res2 ...]
#   Default RES_LIST below.  pixel_size = res * 5 m.  The restart gravitational
#   layer file (hG_ev.asc) is at 35 m, so res*5 >= 35  =>  res >= 7.

set -u

RES_LIST=("${@:-20 15 10 7}")
RES_LIST=(${RES_LIST[@]})

STEPS_SPEED=${STEPS_SPEED:-100}

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

OUTDIR=../Outputs/cpu_gpu_validate
mkdir -p "$OUTDIR"
LOG=$OUTDIR/log.txt
: > "$LOG"

echo "res list: ${RES_LIST[*]}   speed steps: $STEPS_SPEED" | tee -a "$LOG"
echo "date: $(date -Is)" | tee -a "$LOG"
echo | tee -a "$LOG"

# ---------------------------------------------------------------------------
# Pass 1: one-step correctness check
# ---------------------------------------------------------------------------
echo "=== PASS 1: one-step CPU vs GPU correctness check ===" | tee -a "$LOG"
for res in "${RES_LIST[@]}"; do
    echo "--- res=$res (pixel_size=$((res * 5)) m) ---" | tee -a "$LOG"

    STEPS=1 BUILD_DIR=build      bash submission.sh      "$res" 0 >/dev/null 2>&1
    cp "../Outputs/sim_ev_r${res}/run/out.1" "$OUTDIR/cpu1.r${res}.out"

    STEPS=1 BUILD_DIR=build-cuda bash submission_cuda.sh "$res" 0 >/dev/null 2>&1
    cp "../Outputs/sim_ev_r${res}/run/out.1" "$OUTDIR/gpu1.r${res}.out"

    grep '^\[FP\]' "$OUTDIR/cpu1.r${res}.out" > "$OUTDIR/cpu1.r${res}.fp"
    grep '^\[FP\]' "$OUTDIR/gpu1.r${res}.out" > "$OUTDIR/gpu1.r${res}.fp"

    if diff -q "$OUTDIR/cpu1.r${res}.fp" "$OUTDIR/gpu1.r${res}.fp" >/dev/null; then
        echo "  MATCH (byte-identical [FP] lines)" | tee -a "$LOG"
    else
        echo "  MISMATCH:" | tee -a "$LOG"
        diff "$OUTDIR/cpu1.r${res}.fp" "$OUTDIR/gpu1.r${res}.fp" | tee -a "$LOG"
    fi
    paste "$OUTDIR/cpu1.r${res}.fp" "$OUTDIR/gpu1.r${res}.fp" | tee -a "$LOG"
    echo | tee -a "$LOG"
done

# ---------------------------------------------------------------------------
# Pass 2: speedup table
# ---------------------------------------------------------------------------
echo "=== PASS 2: speedup table ($STEPS_SPEED accepted steps) ===" | tee -a "$LOG"
printf "%-8s %-12s %-14s %-14s %-10s\n" "res" "pixel_m" "cpu_wall_s" "gpu_wall_s" "speedup" | tee -a "$LOG"

for res in "${RES_LIST[@]}"; do
    STEPS=$STEPS_SPEED BUILD_DIR=build      bash submission.sh      "$res" 0 >/dev/null 2>&1
    cp "../Outputs/sim_ev_r${res}/run/out.1" "$OUTDIR/cpuN.r${res}.out"

    STEPS=$STEPS_SPEED BUILD_DIR=build-cuda bash submission_cuda.sh "$res" 0 >/dev/null 2>&1
    cp "../Outputs/sim_ev_r${res}/run/out.1" "$OUTDIR/gpuN.r${res}.out"

    cpu_wall=$(grep '^\[BENCH\]' "$OUTDIR/cpuN.r${res}.out" | sed -n 's/.*wall_s=\([0-9.]*\).*/\1/p')
    gpu_wall=$(grep '^\[BENCH\]' "$OUTDIR/gpuN.r${res}.out" | sed -n 's/.*wall_s=\([0-9.]*\).*/\1/p')

    speedup=$(awk -v c="$cpu_wall" -v g="$gpu_wall" 'BEGIN{ if (g>0) printf "%.3f", c/g; else print "n/a" }')

    printf "%-8s %-12s %-14s %-14s %-10s\n" "$res" "$((res * 5))" "$cpu_wall" "$gpu_wall" "$speedup" | tee -a "$LOG"
done

echo | tee -a "$LOG"
echo "done -> $LOG" | tee -a "$LOG"
