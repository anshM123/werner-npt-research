#!/bin/bash
# Runs every check of Sections 3, 5 (Theorems 14', 14'') and 7 of NOTES.md, one Python process at a time.
cd "$(dirname "$0")"
: > logs/run_all.log
for s in selftest check_dilation check_double_encoding check_anchor_persistence check_balanced_projector check_enumerator_lp; do
  start=$(date +%s)
  python $s.py > logs/$s.log 2>&1
  res=$(tail -n 1 logs/$s.log)
  echo "$s : $res ($(( $(date +%s) - start )) s)" | tee -a logs/run_all.log
done
