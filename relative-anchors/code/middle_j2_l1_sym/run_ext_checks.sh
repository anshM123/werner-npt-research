#!/usr/bin/env bash
# Checks for the extension of Theorem L13 to 2 <= b <= 40 (after run_ext30.sh and run_ext40.sh):
# the saved forms of new regimes against the rational-arithmetic reduction, coverage for every k <= 81, and the
# negative controls at sizes inside the new regimes.
cd "$(dirname "$0")"
mkdir -p logs
allok=1
step() {
  name="$1"; exp="$2"; lim="$3"; shift 3
  t0=$(date +%s)
  bash guard.sh "logs/run_${name}.log" "$lim" "$@"
  last=$(tail -n 1 "logs/run_${name}.log" | tr -d '\r')
  if [ "$last" = "$exp" ]; then r=PASS; else r=FAIL; allok=0; fi
  echo "$r  $name  (python $*; expect $exp; $(( $(date +%s) - t0 )) s)"
}
for b in $(seq 14 40); do
  grep -q "^PASS" "logs/ext_cert_${b}.log" && [ "$(tail -n 1 logs/ext_cert_${b}.log)" = PASS ] || { echo "FAIL  ext_cert_${b} not PASS"; allok=0; }
done
step forms_ext PASS 3600 check_forms.py 14+,14 20+,20 27+,27 33+,33 40+,40
step coverage_40_k81 PASS 900 coverage.py 40 none 81
step neg_ext PASS 1800 neg_ext.py
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
