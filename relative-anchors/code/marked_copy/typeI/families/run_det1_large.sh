#!/usr/bin/env bash
# The 19 MT sym families with two symbolic sizes of lists/MT_sym_2e_large.txt (not both i and tau'; the fixed ones
# among i, tau' equal 1; 13 vectors) by Proposition D1 (tdet1.py), spread round-robin over NW worker processes, each
# one Python process under guard.sh (1.2 GB).  Each worker writes its own result pickle res/det1_MT_large_wW.pkl.
# Usage: NW=3 bash run_det1_large.sh   (the logs are logs/run_det1_MT_large_wW.log)
cd "$(dirname "$0")"
NW=${NW:-3}
mkdir -p logs res
fams=($(tr -d '\r' < lists/MT_sym_2e_large.txt))
for w in $(seq 0 $((NW - 1))); do
  mine=()
  for n in "${!fams[@]}"; do
    if [ $((n % NW)) = "$w" ]; then mine+=("${fams[$n]}"); fi
  done
  bash guard.sh logs/run_det1_MT_large_w$w.log 14400 tdet1.py res/det1_MT_large_w$w.pkl "${mine[@]}" &
done
wait
for w in $(seq 0 $((NW - 1))); do echo "worker $w: $(tail -n 1 logs/run_det1_MT_large_w$w.log)"; done
