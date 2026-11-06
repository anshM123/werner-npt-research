#!/usr/bin/env bash
# The saved-block part of the second check at k = 7 (crosscheck.py C) at all 31 values of q from 1e-8 to 1e8, with
# flint evaluation and the exact tests pd_cert / fastpd.py.  Usage: bash run_saved_fast.sh DIR   (DIR holds the cert_k7_*.pkl files)
cd "$(dirname "$0")"
mkdir -p logs
bash guard.sh logs/run_crosscheck_saved_blocks_k7_fast.log 21600 crosscheck.py C "$1" 31
tail -n 1 logs/run_crosscheck_saved_blocks_k7_fast.log
