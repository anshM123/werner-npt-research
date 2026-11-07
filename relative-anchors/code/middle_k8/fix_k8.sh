#!/usr/bin/env bash
# Run after run_k8.sh: repeats the two steps of run_k8.sh whose first run did not count.
#  - coverage_k8 read the log of cert_k8_j4m0_0_0 while that certificate was being rerun with 5.2 GB (fix_big.sh;
#    the rerun passed, logs/fix_big.out), so it saw no PASS line for the configuration (0, 0);
#  - builders_k8 was stopped two seconds after its start, before writing anything (log empty).
# Same step function and limits as run_k8.sh.  Usage: bash fix_k8.sh > logs/run_k8_fix.out
cd "$(dirname "$0")"
all=1
step() {   # step NAME EXPECT LIMIT script args...
  local name="$1" expect="$2" lim="$3"; shift 3
  local log="logs/run_${name}.log"
  local t0=$(date +%s)
  bash guard.sh "$log" "$lim" "$@"
  local st=$?
  local last=$(tail -n 1 "$log" | tr -d '\r')
  local t1=$(date +%s)
  local res=FAIL
  if [ "$st" != 3 ]; then
    if [ "$expect" = PASS ] && [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; fi
    if [ "$expect" = FAIL ] && [ "$last" = FAIL ]; then res=PASS; fi
  fi
  [ "$res" = PASS ] || all=0
  echo "$res  $name  (python $*; expect $expect; $((t1 - t0)) s)"
}
FAMS8="j3m0 j3m1 j3m2 j3m3 j4m0 j4m1 j4m2 j4m3 j4m4 j5m0 j5m1 j5m2 j5m3 j5m4 j5m5 j6m0 j6m1 j6m2 j6m3 j6m4 j6m5 j6m6 j7m2 j7m3 j7m4 j7m5"
ALL8=$(echo $FAMS8 | tr ' ' ',')
echo "fix_k8.sh started $(date)"
grep -q "^PASS$" <(tail -n 1 logs/run_cert_k8_j4m0_0_0.log | tr -d '\r') && echo "cert_k8_j4m0_0_0: PASS (rerun with 5.2 GB, logs/fix_big.out)"
step coverage_k8   PASS 21600 coverage.py "$ALL8" k=8 sample=200
step builders_k8   PASS 21600 check_builders.py "$ALL8" k=8 reps=1
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
