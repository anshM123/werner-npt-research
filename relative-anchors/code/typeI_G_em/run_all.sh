#!/usr/bin/env bash
# Reruns every step of NOTES.md, one Python process at a time, each under guard.sh (memory limit 2 500 000 KB and a
# time limit).  Prints PASS or FAIL per step and ALL PASS at the end.  Usage: bash run_all.sh > logs/run_all.out 2>&1
cd "$(dirname "$0")" || exit 1
mkdir -p logs
echo "run_all.sh started $(date)"
echo "python $(python -c 'import sys; print(sys.version.split()[0])'), python-flint $(python -c 'import flint; print(flint.__version__)'), sympy $(python -c 'import sympy; print(sympy.__version__)'), numpy $(python -c 'import numpy; print(numpy.__version__)')"
allok=1
step () {
  name="$1"; lim="$2"; want="$3"; shift 3
  t0=$(date +%s)
  bash ./guard.sh "logs/run_${name}.log" "$lim" "$@"
  st=$?
  t1=$(date +%s)
  last=$(tail -n 1 "logs/run_${name}.log" | tr -d '\r')
  if [ "$st" = 0 ] && [ "$last" = "$want" ]; then res=PASS; else res=FAIL; allok=0; fi
  echo "step ${name}: ${res} (exit status ${st}, $((t1 - t0)) s, last line: ${last})"
}
step equal_m 1200 PASS equal_m.py
step check_equal_m 1200 PASS check_equal_m.py
step prove_em 600 PASS prove_em.py
step prove_em_H2_5 600 PASS prove_em.py 2 5
step prove_em_H29_100 600 PASS prove_em.py 29 100
step prove_em_H7_10 600 PASS prove_em.py 7 10
step check_e2e 1200 PASS check_e2e.py
step controls_em 600 PASS controls_em.py
step witness_moments 300 PASS witness_moments.py
step ex_frame 900 PASS ex_frame.py
step sweep_general 900 PASS sweep_general.py
step ex_Agrowth 900 "exploration done" ex_Agrowth.py
echo "run_all.sh finished $(date)"
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
