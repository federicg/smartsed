#!/bin/bash
# Nsight Systems profile of the time loop only (cudaProfilerStart/Stop in main). Usage: prof.sh <res> [steps]
source ~/HDD/miniforge3/etc/profile.d/conda.sh; conda activate base
r=$1; n=${2:-100}; P=~/HDD/perf; o=$P/loop_r$r
cd ~/HDD/smartsed/Outputs/sim_ev_r$r/run || exit 1
nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader > $o.gpu_before
nsys profile --capture-range=cudaProfilerApi --capture-range-end=stop --trace=cuda --sample=none --cpuctxsw=none \
  --stats=false --force-overwrite true -o $o $P/main_prof.exe -f SMARTSED_input -sim 0 -scale $r -steps $n > $o.out 2>&1
nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader > $o.gpu_after
for rep in cuda_gpu_kern_sum cuda_api_sum cuda_gpu_mem_time_sum cuda_gpu_mem_size_sum; do
  nsys stats -q --force-export true -r $rep --format csv -o $o $o.nsys-rep >> $o.stats.log 2>&1
done
echo PROF_DONE >> $o.stats.log
