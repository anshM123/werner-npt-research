#!/usr/bin/env bash
# Theorem 7.29.  "./run_checks.sh" verifies the stored certificate; "./run_checks.sh all" first regenerates it.
cd "$(dirname "$0")"
mkdir -p logs cert
if [ "$1" = "all" ]; then
  python tc3_setup.py > logs/tc3_setup.log 2>&1; tail -n 2 logs/tc3_setup.log
  python tc3_round.py > logs/tc3_round.log 2>&1; tail -n 1 logs/tc3_round.log
  python tc3_export.py
fi
python verify_tc3.py > logs/verify_tc3.log 2>&1; tail -n 1 logs/verify_tc3.log
python check_float.py > logs/check_float.log 2>&1; tail -n 1 logs/check_float.log
python tc3_numerics.py > logs/tc3_numerics.log 2>&1; tail -n 1 logs/tc3_numerics.log
python kappa_scan.py 2 restricted none sn CLARABEL > logs/kappa_scan_n2.log 2>&1; tail -n 1 logs/kappa_scan_n2.log
python kappa_scan.py 3 restricted none sn CLARABEL > logs/kappa_scan_n3.log 2>&1; tail -n 1 logs/kappa_scan_n3.log
