#!/usr/bin/env bash
# The second, independent checks of crosscheck.py (one Python process under guard.sh; PASS or FAIL, ALL PASS at the end).
# Usage: bash run_crosscheck.sh > logs/run_crosscheck.out
cd "$(dirname "$0")"
mkdir -p logs
t0=$(date +%s)
echo "run_crosscheck.sh started $(date)"
bash guard.sh logs/run_crosscheck.log 3600 crosscheck.py all
st=$?
last=$(tail -n 1 logs/run_crosscheck.log | tr -d '\r')
if [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; else res=FAIL; fi
echo "$res  crosscheck  (python crosscheck.py all; $(( $(date +%s) - t0 )) s)"
echo "run_crosscheck.sh finished $(date)"
if [ "$res" = PASS ]; then echo "ALL PASS"; else echo "NOT ALL PASS"; fi
