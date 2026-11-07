#!/usr/bin/env bash
# Every check of README.md; logs in ../logs.  Needs Node (no packages) and Python with NumPy and SciPy.
cd "$(dirname "$0")"
mkdir -p ../logs
all=1
python run_expr.py transfer_certificates.js > ../logs/transfer_certificates.log 2>&1 && grep -q '"passed": true' ../logs/transfer_certificates.log && echo "transfer certificates: PASS" || { echo "transfer certificates: FAIL"; all=0; }
python run_expr.py tail_certificates.js > ../logs/tail_certificates.log 2>&1 && grep -q '"passed": true' ../logs/tail_certificates.log && echo "tail certificates: PASS" || { echo "tail certificates: FAIL"; all=0; }
python cyclic_check.py > ../logs/cyclic_check.log 2>&1 && cat ../logs/cyclic_check.log || { echo "cyclic_check: FAIL"; all=0; }
# sizes 14..17 and 27..29 (about 15 minutes)
gap() { local log="../logs/gap_$1.log"; shift; python gap_certificates.py "$@" > "$log" 2>&1; local v=$(tail -n 1 "$log" | tr -d '\r'); echo "gap_certificates $*: $v"; [ "$v" = PASS ] || all=0; }
gap control control
gap rows_12 rows 12
for nl in "14 7" "15 8" "16 8" "17 9"; do set -- $nl; gap split_$1 split $1 $2; done
for L2 in 14 15 16; do gap block_13_$L2 block 13 $L2; done
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
