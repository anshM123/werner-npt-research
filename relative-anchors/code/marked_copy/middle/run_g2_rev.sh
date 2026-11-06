#!/usr/bin/env bash
# Extra workers for the same grid, taking the columns from the top (y = 98 downwards); slices already saved are
# skipped, so the workers of run_g2.sh and these share the grid.
cd "$(dirname "$0")"
NW=${NW:-3}
for w in $(seq 0 $((NW - 1))); do
  G2_REVERSE=1 bash guard.sh logs/g2_rworker_$w.log 36000 g2slices.py save g2out 0 98 0 102 $w $NW &
done
wait
for w in $(seq 0 $((NW - 1))); do echo "reverse worker $w: $(tail -n 1 logs/g2_rworker_$w.log)"; done
