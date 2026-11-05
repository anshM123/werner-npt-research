#!/bin/bash
# Runs every check of NOTES.md one after the other (one Python process at a time).
cd "$(dirname "$0")"
: > run_all.log
for s in check_closed_form check_exact sweep; do
  start=$(date +%s)
  python $s.py > $s.log 2>&1
  res=$(tail -n 1 $s.log)
  echo "$s : $res ($(( $(date +%s) - start )) s)" | tee -a run_all.log
done
