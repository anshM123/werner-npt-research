#!/usr/bin/env bash
# Every check of this folder, with logs in logs/.  Needs Node (no packages) and Python with NumPy.
cd "$(dirname "$0")"
python run_js.py diag_certificate.js sfP6Verify > logs/diag_certificate.log 2>&1; echo "diag_certificate: $(head -c 400 logs/diag_certificate.log | grep -o '"status": "[A-Z]*"')"
python control_H.py > logs/control_H.log 2>&1; echo "control_H: $(tail -n 1 logs/control_H.log)"
python diag_check.py > logs/diag_check.log 2>&1; tail -n 1 logs/diag_check.log
python offdiag_ratio.py > logs/offdiag_ratio.log 2>&1; head -n 1 logs/offdiag_ratio.log
rm -rf __pycache__
