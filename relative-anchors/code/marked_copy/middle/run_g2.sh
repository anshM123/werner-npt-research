#!/usr/bin/env bash
# The slices of the determinant route on the full verification grid: y = 0..98 (t = 4..102), z = 0..102
# (i = 7..109), NW workers in parallel (each one Python process under guard.sh, 2.5 GB), then stage B.
cd "$(dirname "$0")"
NW=${NW:-5}
mkdir -p g2out logs
for w in $(seq 0 $((NW - 1))); do
  bash guard.sh logs/g2_worker_$w.log 36000 g2slices.py save g2out 0 98 0 102 $w $NW &
done
wait
for w in $(seq 0 $((NW - 1))); do echo "worker $w: $(tail -n 1 logs/g2_worker_$w.log)"; done
bash guard.sh logs/g2_stageB.log 7200 g2interp.py B g2out 98 56 60 200
echo "stage B: $(tail -n 1 logs/g2_stageB.log)"
