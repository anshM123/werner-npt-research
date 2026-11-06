#!/usr/bin/env bash
# The end-to-end check at k = 7 with the faster exact routines (NOTES.md Section 3), one Python process at a time.
# Usage: bash run_end2end_fast.sh   (writes logs/validate_*.log and logs/run_end2end_k7_fast.log)
cd "$(dirname "$0")"
mkdir -p logs
bash guard.sh logs/validate_pd.log 1800 validate_pd.py 30 11
bash guard.sh logs/validate_fastj1_k4.log 1800 validate.py 4 2 2
bash guard.sh logs/validate_fastj1_k5.log 3600 validate.py 5 2 3
bash guard.sh logs/validate_fastj1_k6.log 7200 validate.py 6 1 4
bash guard.sh logs/run_end2end_k7_fast.log 21600 check_end2end.py j3m0,j3m1,j3m2,j3m3,j4m0,j4m1,j4m2,j4m3,j4m4,j5m0,j5m1,j5m2,j5m3,j5m4,j5m5,j6m2,j6m3,j6m4 k=7 reps=1
for f in logs/validate_pd.log logs/validate_fastj1_k4.log logs/validate_fastj1_k5.log logs/validate_fastj1_k6.log logs/run_end2end_k7_fast.log; do echo "$(tail -n 1 $f | tr -d '\r')  $f"; done
