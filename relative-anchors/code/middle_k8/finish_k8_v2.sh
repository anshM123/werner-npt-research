#!/usr/bin/env bash
# The saved-blocks step was started twice under guard_big.sh. The second start was stopped by hand after 2.5 hours on
# cert_k8_j4m0_0_4 (its 40 earlier blocks had passed; logs/stopped_run_crosscheck_saved.log): at q = 10^8 the
# floating Cholesky step of pd_cert fails on that block, and the exact LDL fallback is very slow at that size.
# crosscheck.py now tries pd_cert2 (the same exact test with a sturdier floating step) before the exact LDL.
# This script reruns, one after the other under guard_big.sh: the negative controls (now also for pd_cert2), the
# whole saved-blocks step, and the nine end-to-end families that the 2.5 GB limit had stopped.
# Usage: bash finish_k8_v2.sh > logs/run_k8_v2.out
cd "$(dirname "$0")"
step() {   # step NAME LOG SECONDS args...
  name="$1"; lg="$2"; secs="$3"; shift 3
  t0=$(date +%s)
  bash guard_big.sh "$lg" "$secs" "$@"
  st=$?; last=$(tail -n 1 "$lg" | tr -d '\r')
  if [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; else res=FAIL; fi
  echo "$res  $name  (python $*, under guard_big.sh; expect PASS; $(( $(date +%s) - t0 )) s)"
}
step crosscheck_neg logs/run_crosscheck_neg.log 3600 crosscheck.py N | tee logs/run_k8_neg.out
step crosscheck_saved logs/run_crosscheck_saved.log 43200 crosscheck.py C logs 31 | tee logs/run_k8_saved.out
REST="j4m0,j4m1,j4m2,j4m3,j4m4,j5m0,j5m1,j5m2,j5m3"
step end2end_first_k8_h1rest logs/run_end2end_first_k8_h1rest.log 43200 check_end2end.py "$REST" k=8 reps=1 first=1 | tee logs/run_k8_h1rest.out
