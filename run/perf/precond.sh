#!/bin/bash
# IC(0) vs Jacobi vs plain CG: time loop only ([BENCH] timer), 100 accepted steps. Usage: precond.sh <res>...
source ~/HDD/miniforge3/etc/profile.d/conda.sh; conda activate base
P=~/HDD/perf
for r in "$@"; do
  cd ~/HDD/smartsed/Outputs/sim_ev_r$r/run || exit 1
  for v in ic0 jacobi none; do
    case $v in
      ic0)    sed "s/^\( *use_preconditioner *=\) *[a-z]*.*/\1 true\n    jacobi = false/" SMARTSED_input > in_$v ;;
      jacobi) sed "s/^\( *use_preconditioner *=\) *[a-z]*.*/\1 true\n    jacobi = true/" SMARTSED_input > in_$v ;;
      none)   sed "s/^\( *use_preconditioner *=\) *[a-z]*/\1 false/" SMARTSED_input > in_$v ;;
    esac
    out=$P/pc_r${r}_$v.out
    echo "gpu_util_before=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader)" > $out
    $P/main_jac.exe -f in_$v -sim 0 -scale $r -steps 100 >> $out 2>&1
    echo "res=$r $v $(grep -o "loop_iters=[0-9]*  wall_s=[0-9.]*.*" $out) cg_iters=$(grep -o "iterations performed: [0-9]*" $out | awk "{s+=\$3} END{print s}") cg_max=$(grep -o "iterations performed: [0-9]*" $out | awk "{if(\$3>m)m=\$3} END{print m}") $(grep "^\[FP\] H " $out | tail -1) $(head -1 $out)"
  done
done
echo PC_DONE
