#!/usr/bin/env bash
# The regime of Theorem F1b left open by run_all.sh and run_fext.sh: family 1 (T = {t, tau}) with t and tau of different
# weights (tw = (0, 1)) and both U-levels symbolic (2+, 2+): the invariant form (8 leading minors) and the two standard
# forms; then the coverage with this regime included.  One Python process at a time under guard.sh (2.5 GB).
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
echo "run_open.sh started $(date)"
step f1_open_std0 PASS 3600 cert_t2.py 2+ 2+ part=std0 tw=0,1
step f1_open_std1 PASS 3600 cert_t2.py 2+ 2+ part=std1 tw=0,1
step f1_open_inv PASS 21600 cert_t2.py 2+ 2+ part=inv tw=0,1
step coverage_open PASS 900 coverage.py 40 M1=40 KMAX=120 F1OPEN=1
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
