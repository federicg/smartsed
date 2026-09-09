#!/bin/bash
#
# Full time-loop speedup sweep (buildMatrix + IC(0)-preconditioned PCG solve
# active, red-black reordering on the GPU side, updateVel/sediment still
# disabled) across all 8 resolutions. One accepted step each (-steps 1),
# CPU vs GPU, both backends' wall-clock time and H/eta fingerprints.

set -u

RES_LIST=("${@:-20 15 10 7 5 3 2 1}")
RES_LIST=(${RES_LIST[@]})

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

OUTDIR=../Outputs/speedup_full_solve
mkdir -p "$OUTDIR"
LOG=$OUTDIR/log.txt
: > "$LOG"

echo "res list: ${RES_LIST[*]}   (full time loop: buildMatrix + PCG solve, -steps 1)" | tee -a "$LOG"
echo "date: $(date -Is)" | tee -a "$LOG"
echo | tee -a "$LOG"

printf "%-6s %-10s %-12s %-14s %-14s %-10s %-14s %-14s\n" \
  "res" "pixel_m" "basin_cells" "cpu_wall_s" "gpu_wall_s" "speedup" "cpu_H_sum" "gpu_H_sum" | tee -a "$LOG"

for res in "${RES_LIST[@]}"; do
    STEPS=1 BUILD_DIR=build      bash submission.sh      "$res" 0 >/dev/null 2>&1
    cp "../Outputs/sim_ev_r${res}/run/out.1" "$OUTDIR/cpu.r${res}.out"

    STEPS=1 BUILD_DIR=build-cuda bash submission_cuda.sh "$res" 0 >/dev/null 2>&1
    cp "../Outputs/sim_ev_r${res}/run/out.1" "$OUTDIR/gpu.r${res}.out"

    cpu_wall=$(grep '^\[BENCH\]' "$OUTDIR/cpu.r${res}.out" | sed -n 's/.*wall_s=\([0-9.]*\).*/\1/p')
    gpu_wall=$(grep '^\[BENCH\]' "$OUTDIR/gpu.r${res}.out" | sed -n 's/.*wall_s=\([0-9.]*\).*/\1/p')
    ncells=$(grep '^\[FP\] H' "$OUTDIR/cpu.r${res}.out" | sed -n 's/.*n=\([0-9]*\).*/\1/p')
    cpu_H=$(grep '^\[FP\] H' "$OUTDIR/cpu.r${res}.out" | sed -n 's/.*sum=\([+-][0-9.eE+-]*\).*/\1/p')
    gpu_H=$(grep '^\[FP\] H' "$OUTDIR/gpu.r${res}.out" | sed -n 's/.*sum=\([+-][0-9.eE+-]*\).*/\1/p')

    speedup=$(awk -v c="$cpu_wall" -v g="$gpu_wall" 'BEGIN{ if (g>0) printf "%.3f", c/g; else print "n/a" }')

    printf "%-6s %-10s %-12s %-14s %-14s %-10s %-14s %-14s\n" \
      "$res" "$((res * 5))" "$ncells" "$cpu_wall" "$gpu_wall" "$speedup" "$cpu_H" "$gpu_H" | tee -a "$LOG"
done

echo | tee -a "$LOG"
echo "done -> $LOG" | tee -a "$LOG"
