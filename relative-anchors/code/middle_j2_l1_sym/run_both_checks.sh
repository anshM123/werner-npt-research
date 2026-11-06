#!/usr/bin/env bash
# Checks for the regime (14+, 14+) (both level sizes symbolic): every family certified (logs/both_*.log end in PASS),
# the saved forms against the rational-arithmetic reduction, and the coverage of every two-value configuration with
# 11 <= k <= 120 by the regimes (o(m)+, m), m <= 13, and (14+, 14+).
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
for f in inv mix std0 two0 anti0 three0 tto0 twostd0 antistd0 std1 two1 anti1 three1 tto1 twostd1 antistd1; do
  if [ "$(tail -n 1 logs/both_$f.log | tr -d '\r')" = PASS ]; then echo "PASS  both_$f  ($(grep '^time' logs/both_$f.log))"; else echo "FAIL  both_$f"; allok=0; fi
done
step forms_collect PASS 600 build_forms_both.py
step forms_both PASS 7200 check_forms.py 14+,14+
step coverage_both_k120 PASS 1800 coverage.py 13 14 120
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
