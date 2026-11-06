#!/usr/bin/env bash
# Rerun of the builder step of run_k7.sh with the corrected check_builders.py (NOTES.md Section 4).
# Usage: bash run_k7_builders.sh > logs/run_k7_builders.out   (needs the blocks saved by run_k7.sh in logs/)
cd "$(dirname "$0")"
mkdir -p logs
FAMS7="j3m0 j3m1 j3m2 j3m3 j4m0 j4m1 j4m2 j4m3 j4m4 j5m0 j5m1 j5m2 j5m3 j5m4 j5m5 j6m2 j6m3 j6m4"
ALL7=$(echo $FAMS7 | tr ' ' ',')
t0=$(date +%s)
echo "run_k7_builders.sh started $(date)"
log="logs/run_k7_builders_builders_k7.log"
bash guard.sh "$log" 21600 check_builders.py "$ALL7" k=7 reps=1
st=$?
last=$(tail -n 1 "$log" | tr -d '\r')
res=FAIL
if [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; fi
t1=$(date +%s)
echo "$res  builders_k7  (python check_builders.py $ALL7 k=7 reps=1; expect PASS; $((t1 - t0)) s)"
if [ "$res" = PASS ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
