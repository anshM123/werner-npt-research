#!/usr/bin/env bash
# The second extension (Section 8 of NOTES.md), run after run_extra.sh: the fixed-k certificate for k = 29 and
# the coverage and block bookkeeping for 8 <= k <= 29 (with the files of run_all.sh and run_extra.sh).  One Python
# process at a time under guard.sh; PASS or FAIL per step, ALL PASS at the end.
# Usage: bash run_extra2.sh > logs/run_extra2.out
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

echo "run_extra2.sh started $(date)"
for k in 29; do
  step cert_$k        PASS 3600 cert_t.py $k $k 1 logs/cert_$k.pkl
done
PK="logs/cert_08_12.pkl logs/cert_13_15.pkl logs/cert_16_17.pkl logs/cert_18.pkl logs/cert_19.pkl logs/cert_20.pkl"
for k in 21 22 23 24 25 26 27 28 29; do PK="$PK logs/cert_$k.pkl"; done
step coverage_8_29    PASS 1800 coverage_t.py 8 29 $PK
step blocks_8_29      PASS 1800 blocks_t.py 8 8 29 $PK
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
