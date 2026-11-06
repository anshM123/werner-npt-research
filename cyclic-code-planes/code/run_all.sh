#!/usr/bin/env bash
# Every check of README.md; logs in ../logs.  Needs Node (no packages) and Python with NumPy and SciPy.
cd "$(dirname "$0")"
mkdir -p ../logs
all=1
python run_expr.py transfer_certificates.js > ../logs/transfer_certificates.log 2>&1 && grep -q '"passed": true' ../logs/transfer_certificates.log && echo "transfer certificates: PASS" || { echo "transfer certificates: FAIL"; all=0; }
python run_expr.py tail_certificates.js > ../logs/tail_certificates.log 2>&1 && grep -q '"passed": true' ../logs/tail_certificates.log && echo "tail certificates: PASS" || { echo "tail certificates: FAIL"; all=0; }
python cyclic_check.py > ../logs/cyclic_check.log 2>&1 && cat ../logs/cyclic_check.log || { echo "cyclic_check: FAIL"; all=0; }
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
