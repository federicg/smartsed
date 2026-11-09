#!/bin/bash
# A/B: cublasDnrm2 (old) vs dot+sqrt (new), time loop only ([BENCH] timer), 100 steps.
source ~/HDD/miniforge3/etc/profile.d/conda.sh; conda activate base
while pgrep -u fgatti -f "nsys|main.exe" >/dev/null; do sleep 10; done
for r in 3 2 1; do
  cd ~/HDD/smartsed/Outputs/sim_ev_r$r/run || exit 1
  for v in nrm2 dot; do
    cp ~/HDD/perf/main_$v.exe ../DeterministicProgram/main.exe
    for k in 1 2 3; do
      out=~/HDD/perf/ab_r${r}_${v}_$k.out
      ../DeterministicProgram/main.exe -f SMARTSED_input -sim 0 -scale $r -steps 100 > $out 2>&1
      ms=$(sed -n "s/.*per_iter_ms=\([0-9.]*\).*/\1/p" $out)
      cg=$(grep -o "iterations performed: [0-9]*" $out | awk "{s+=\$3} END{print s}")
      fp=$(grep "^\[FP\] H" $out | tail -1 | sed -n "s/.*sum=\([^ ]*\).*/\1/p")
      echo "res=$r $v run=$k per_iter_ms=$ms cg_iters=$cg H_sum=$fp"
    done
  done
done
echo AB_DONE
