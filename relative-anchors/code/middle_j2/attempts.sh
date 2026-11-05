#!/usr/bin/env bash
# The two regimes whose certificates did not finish within the limits (NOTES.md, rows 5b and 10b).  Not part of
# run_all.sh.  In the recorded attempt the first one passed 7 of its 9 leading minors and was stopped by the 2.5 GB
# limit during minor 8; the second passed 7 of its 8 leading minors and was stopped by hand during minor 8
# (logs/attempts/stopped_regimes.txt).
# Usage: bash attempts.sh
cd "$(dirname "$0")"
mkdir -p logs/attempts
bash guard.sh logs/attempts/f2_inv_3_3.log 3600 cert_f2.py inv 3+ 3+ tw=0
echo "f2_inv_3_3: exit $?, last line: $(tail -n 1 logs/attempts/f2_inv_3_3.log)"
bash guard.sh logs/attempts/f1_cert_B_2_2.log 5400 cert_t2.py 2+ 2+ tw=0,1
echo "f1_cert_B_2_2: exit $?, last line: $(tail -n 1 logs/attempts/f1_cert_B_2_2.log)"
