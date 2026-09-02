#!/bin/bash
#
# CPU-vs-GPU speedup sweep for the SMARTSED CUDA port.
#
# For each DEM coarsening factor in RES_LIST it runs the CPU-Release build and
# the CUDA build, each stopping after $STEPS accepted time steps (passed as
# "-steps N" through the submission scripts) and collects the [BENCH]
# wall-clock line and the [FP] field-fingerprint lines that main_final_H.cpp
# prints at the end of the time loop.
#
# Usage:  ./benchmark_sweep.sh [res1 res2 ...]
#   default RES_LIST below.  pixel_size = res * 5 m.
#   Env: STEPS (default 100).  The GPU CG solve currently aborts at res=2
#        (10 m), so it is kept out of the default list.
#
# Results:  ../Outputs/bench_sweep/summary.txt  (+ per-run out/err kept in
#           ../Outputs/sim_ev_r<res>/run for the last run of each resolution)

set -u

RES_LIST=("${@:-20 10 5 3}")
# allow a single quoted "20 10 5 3" as well as separate args
RES_LIST=(${RES_LIST[@]})

STEPS=${STEPS:-100}                    # accepted time steps per run
export STEPS

CPU_BUILD=${CPU_BUILD:-build}          # CPU Release
GPU_BUILD=${GPU_BUILD:-build-cuda}     # CUDA

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

OUTDIR=../Outputs/bench_sweep
mkdir -p "$OUTDIR"
SUMMARY=$OUTDIR/summary.txt
: > "$SUMMARY"

echo "res list: ${RES_LIST[*]}   steps: $STEPS"        | tee -a "$SUMMARY"
echo "cpu build: $CPU_BUILD   gpu build: $GPU_BUILD"   | tee -a "$SUMMARY"
echo "date: $(date -Is)"                               | tee -a "$SUMMARY"
echo                                                   | tee -a "$SUMMARY"

run_one () {           # $1 = tag (CPU|GPU)  $2 = script  $3 = build dir  $4 = res
    local tag=$1 script=$2 build=$3 res=$4
    echo "=== $tag  res=$res  (pixel_size=$((res * 5)) m) ===" | tee -a "$SUMMARY"
    BUILD_DIR=$build bash "$script" "$res" 0
    local out=../Outputs/sim_ev_r${res}/run/out.1
    cp "$out" "$OUTDIR/out.${tag}.r${res}" 2>/dev/null
    grep -E '^\[BENCH\]|^\[FP\]|cell resolution|dt_DSV_given' "$out" \
        | sed "s/^/  $tag r$res | /" | tee -a "$SUMMARY"
    echo | tee -a "$SUMMARY"
}

for res in "${RES_LIST[@]}"; do
    run_one CPU submission.sh      "$CPU_BUILD" "$res"
    run_one GPU submission_cuda.sh "$GPU_BUILD" "$res"
done

echo "done -> $SUMMARY"
