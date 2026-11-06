#!/usr/bin/env bash
# The second check (crosscheck.py) at k = 7, one Python process at a time under guard.sh.
# Usage: bash run_crosscheck_k7.sh > logs/run_crosscheck_k7.out   (the last step reads the blocks saved by run_k7.sh)
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
echo "run_crosscheck_k7.sh started $(date)"
step negative_controls  900   N
step all_blocks_k7      14400 A 7 24 7207
step saved_blocks_k7    14400 C logs 16
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
