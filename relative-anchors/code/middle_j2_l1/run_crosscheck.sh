#!/usr/bin/env bash
# The second check of Section 8 of NOTES.md (crosscheck.py), one Python process at a time under guard.sh.
# Usage: bash run_crosscheck.sh > logs/run_crosscheck.out   (mode C reads the forms saved by run_all.sh in logs/)
cd "$(dirname "$0")"
mkdir -p logs
all=1
t_start=$(date +%s)
step() {   # step NAME LIMIT args...
  local name="$1" lim="$2"; shift 2
  local log="logs/run_crosscheck_${name}.log"
  local t0=$(date +%s)
  bash guard.sh "$log" "$lim" crosscheck.py "$@"
  local st=$?
  local last=$(tail -n 1 "$log" | tr -d '\r')
  local t1=$(date +%s)
  local res=FAIL
  if [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; fi
  [ "$res" = PASS ] || all=0
  echo "$res  $name  (python crosscheck.py $*; $((t1 - t0)) s)"
}
echo "run_crosscheck.sh started $(date)"
step negative_controls 1800 N
step ylib              1800 B
step saved_forms       1800 C logs
step end_to_end_k5_k8  3600 A 5 8
step end_to_end_k9     3600 A 9 9
step end_to_end_k10    7200 A 10 10
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
