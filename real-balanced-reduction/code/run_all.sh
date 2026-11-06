#!/usr/bin/env bash
# Every check of README.md; logs in ../logs.  Needs Node (no packages) and Python with NumPy.
cd "$(dirname "$0")"
mkdir -p ../logs
all=1
python run_js.py realification_exact.js verifyRankPreservingRealification > ../logs/realification_exact.log 2>&1 && echo "realification (exact): PASS" || { echo "realification (exact): FAIL"; all=0; }
python run_js.py balancing_exact.js verifySymmetryInterpolation > ../logs/balancing_exact.log 2>&1 && echo "balancing (exact): PASS" || { echo "balancing (exact): FAIL"; all=0; }
python realify_check.py > ../logs/realify_check.log 2>&1; tail -n 1 ../logs/realify_check.log | grep -q PASS && echo "realification identity (from the definition): PASS" || { echo "realify_check: FAIL"; all=0; }
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
