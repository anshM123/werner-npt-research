#!/usr/bin/env bash
# Reruns every check of NOTES.md, one Python process at a time, each under guard.sh (memory limit 2 500 000 KB and a
# time limit).  Prints PASS or FAIL per step and ALL PASS at the end.  Usage: bash run_all.sh > logs/run_all.out 2>&1
cd "$(dirname "$0")" || exit 1
mkdir -p logs
echo "run_all.sh started $(date)"
echo "python $(python -c 'import sys; print(sys.version.split()[0])'), python-flint $(python -c 'import flint; print(flint.__version__)'), numpy $(python -c 'import numpy; print(numpy.__version__)')"
allok=1
step () {
  name="$1"; lim="$2"; want="$3"
  t0=$(date +%s)
  ./guard.sh "logs/${name}.log" "$lim" "${name}.py"
  st=$?
  t1=$(date +%s)
  last=$(tail -n 1 "logs/${name}.log" | tr -d '\r')
  if [ "$st" = 0 ] && [ "$last" = "$want" ]; then res=PASS; else res=FAIL; allok=0; fi
  echo "step ${name}: ${res} (exit status ${st}, $((t1 - t0)) s, last line: ${last})"
}
step check_closed 1800 PASS
step cert_pair 600 PASS
step cert_row 900 PASS
step check_reduction 1800 PASS
step check_pd 900 PASS
step negative_control 600 PASS
step independent_cert 1200 PASS
step independent_e2e 1800 PASS
step numerics 900 "numerical illustration done"
echo "run_all.sh finished $(date)"
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
