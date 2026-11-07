#!/usr/bin/env bash
# Two steps of finish_k8.sh were stopped by the 2.5 GB limit of guard.sh: crosscheck_saved (at the 616-state blocks of
# j4m0; every block before them had passed) and end2end_first_k8_h1 (after its first four families j3m0..j3m3 had
# passed).  This script reruns, one after the other under guard_big.sh (5.2 GB), the whole saved-blocks step and the
# nine remaining end-to-end families.  Usage: bash finish_k8_rest.sh > logs/run_k8_rest.out
cd "$(dirname "$0")"
t0=$(date +%s)
bash guard_big.sh logs/run_crosscheck_saved.log 43200 crosscheck.py C logs 31
st=$?; last=$(tail -n 1 logs/run_crosscheck_saved.log | tr -d '\r')
if [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; else res=FAIL; fi
echo "$res  crosscheck_saved  (python crosscheck.py C logs 31, under guard_big.sh; expect PASS; $(( $(date +%s) - t0 )) s)" | tee logs/run_k8_saved.out
REST="j4m0,j4m1,j4m2,j4m3,j4m4,j5m0,j5m1,j5m2,j5m3"
t0=$(date +%s)
bash guard_big.sh logs/run_end2end_first_k8_h1rest.log 43200 check_end2end.py "$REST" k=8 reps=1 first=1
st=$?; last=$(tail -n 1 logs/run_end2end_first_k8_h1rest.log | tr -d '\r')
if [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; else res=FAIL; fi
echo "$res  end2end_first_k8_h1rest  (python check_end2end.py $REST k=8 reps=1 first=1, under guard_big.sh; expect PASS; $(( $(date +%s) - t0 )) s)" | tee logs/run_k8_h1rest.out
