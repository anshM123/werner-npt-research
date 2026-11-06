#!/usr/bin/env bash
# The reduced forms against the Lemma B block at k = 9 (every j) and k = 10 (j = 3, 8, 9), where all orbit types
# occur together (m >= 2, i >= 2 and tau' >= 2 at once), run after run_extra.sh.  One Python process at a time under
# guard.sh; PASS or FAIL per step, ALL PASS at the end.
# Usage: bash run_check2.sh > logs/run_check2.out
cd "$(dirname "$0")"
mkdir -p logs
all=1
t_start=$(date +%s)

step() {   # step NAME EXPECT LIMIT script args...
  local name="$1" expect="$2" lim="$3"; shift 3
  local log="logs/run_${name}.log"
  local t0=$(date +%s)
  bash guard.sh "$log" "$lim" "$@"
  local st=$?
  local last=$(tail -n 1 "$log" | tr -d '\r' | awk '{print $1}')
  local t1=$(date +%s)
  local res=FAIL
  if [ "$st" != 3 ]; then
    if [ "$expect" = PASS ] && [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; fi
    if [ "$expect" = FAIL ] && [ "$last" = FAIL ]; then res=PASS; fi
  fi
  [ "$res" = PASS ] || all=0
  echo "$res  $name  (python $*; expect $expect; $((t1 - t0)) s)"
}

echo "run_check2.sh started $(date)"
step explicit_9       PASS 3000 check_explicit.py 9 9 1 3/2
step explicit_10      PASS 3000 check_explicit.py 10 10 3 1/2 3,8,9
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
