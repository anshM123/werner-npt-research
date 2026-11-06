#!/usr/bin/env bash
# The second check for the parts of irrep index i >= 5 (middle/NOTES_MI.md): mfcheck.py at 60 members per region with
# 56 <= k <= 100, beyond every fixed-k certificate.  One Python process at a time under guard.sh.
# Usage: bash run_check_mi.sh > logs/run_check_mi.out
cd "$(dirname "$0")"
mkdir -p logs
all=1
step() {   # step NAME LIMIT script args...
  local name="$1" lim="$2"; shift 2
  local log="logs/run_check_${name}.log"
  local t0=$(date +%s)
  bash guard.sh "$log" "$lim" "$@"
  local st=$?
  local last=$(tail -n 1 "$log" | tr -d '\r')
  local res=FAIL
  if [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; fi
  [ "$res" = PASS ] || all=0
  echo "$res  $name  (python $*; $(( $(date +%s) - t0 )) s)"
}
echo "run_check_mi.sh started $(date)"
step members_mi_n56  1800 mfcheck.py n56 60 501 100
step members_mi_eq   1800 mfcheck.py eq 60 502 100
step members_mi_et   1800 mfcheck.py et 60 503 100
step members_mi_g7   1800 mfcheck.py g7 60 504 100
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
