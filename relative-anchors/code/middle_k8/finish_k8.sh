#!/usr/bin/env bash
# Run after run_k8.sh, which was stopped during its step neg_points_k8: at k = 8 that step needs about 18 hours (32 of
# its 240 configurations in 2.4 hours) against a limit of 6, and check_end2end.py would need about 30.  The log of the
# stopped step is logs/stopped_run_neg_points_k8.log.  This script runs what was left, with the same step function:
#  lane 1: the second check (crosscheck.py) in full -- the negative controls, every torus block of every (R_inf)_j at
#          12 random two-value points (three runs of 4 points, seeds 8208, 8209, 8210), the saved blocks at 31 q;
#  lane 2 (after the other long job on this machine has ended): neg_points.py and check_end2end.py on the first
#          two-level representative configuration of every family (first=1), the end-to-end check in two halves.
# Usage: bash finish_k8.sh > logs/run_k8_finish.out
cd "$(dirname "$0")"
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
  echo "$res  $name  (python $*; expect $expect; $((t1 - t0)) s)"
}
FAMS8="j3m0 j3m1 j3m2 j3m3 j4m0 j4m1 j4m2 j4m3 j4m4 j5m0 j5m1 j5m2 j5m3 j5m4 j5m5 j6m0 j6m1 j6m2 j6m3 j6m4 j6m5 j6m6 j7m2 j7m3 j7m4 j7m5"
ALL8=$(echo $FAMS8 | tr ' ' ',')
HALF1="j3m0,j3m1,j3m2,j3m3,j4m0,j4m1,j4m2,j4m3,j4m4,j5m0,j5m1,j5m2,j5m3"
HALF2="j5m4,j5m5,j6m0,j6m1,j6m2,j6m3,j6m4,j6m5,j6m6,j7m2,j7m3,j7m4,j7m5"
echo "finish_k8.sh started $(date)"
( step crosscheck_neg       PASS 900   crosscheck.py N
  step crosscheck_blocks_a  PASS 21600 crosscheck.py A 8 4 8208
  step crosscheck_blocks_b  PASS 21600 crosscheck.py A 8 4 8209
  step crosscheck_blocks_c  PASS 21600 crosscheck.py A 8 4 8210
  step crosscheck_saved     PASS 21600 crosscheck.py C logs 31 ) > logs/finish_lane1.out 2>&1 &
( until grep -q "^run_verify_st4.sh done" ../verifyST4/logs/verify_st4.out 2>/dev/null; do sleep 120; done
  step neg_points_first_k8   PASS 21600 neg_points.py "$ALL8" k=8 reps=1 first=1
  step end2end_first_k8_h1   PASS 43200 check_end2end.py "$HALF1" k=8 reps=1 first=1
  step end2end_first_k8_h2   PASS 43200 check_end2end.py "$HALF2" k=8 reps=1 first=1 ) > logs/finish_lane2.out 2>&1 &
wait
cat logs/finish_lane1.out logs/finish_lane2.out
n=$(cat logs/finish_lane1.out logs/finish_lane2.out | grep -c "^PASS ")
if [ "$n" = 8 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED ($n of 8 passed)"; fi
