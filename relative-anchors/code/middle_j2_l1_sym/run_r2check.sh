#!/usr/bin/env bash
# The second check: the whole letter-1 block of (R_inf)_2 from Lemma B at random two-value points, 11 <= k <= 14.
cd "$(dirname "$0")"
mkdir -p logs
for k in 11 12 13; do bash guard.sh logs/r2check_k$k.log 3000 r2check.py $k 3 $((k * 7 + 1)); tail -n 1 logs/r2check_k$k.log; done
bash guard.sh logs/r2check_k14.log 3000 r2check.py 14 2 99; tail -n 1 logs/r2check_k14.log
