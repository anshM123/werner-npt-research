#!/usr/bin/env bash
# Theorems F1b and F2b extended: the fixed-size regimes of run_all.sh for 13 <= m <= 40 (family 2: inv part at (m, 3+)
# and (3+, m); family 1: tw = (0, 1) at (m, 2+)), then coverage with M2 = M1 = 40.  One Python process at a time.
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
for m in $(seq 13 40); do
  step f2_inv_${m}_3p PASS 1800 cert_f2.py inv "$m" 3+ tw=0
  step f2_inv_3p_${m} PASS 1800 cert_f2.py inv 3+ "$m" tw=0
  step f1_cert_B_${m}_2 PASS 1800 cert_t2.py "$m" 2+ tw=0,1
done
step coverage_40 PASS 900 coverage.py 40 M1=40 KMAX=100
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
