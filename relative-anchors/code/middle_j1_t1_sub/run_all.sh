#!/usr/bin/env bash
# Reruns every check of NOTES.md, one Python process at a time, each under guard.sh (memory limit 2 500 000 KB and a
# time limit).  Prints PASS or FAIL per step and ALL PASS at the end.  Usage: bash run_all.sh > logs/run_all.out 2>&1
cd "$(dirname "$0")" || exit 1
export PYTHONWARNINGS=ignore    # library warnings (solver accuracy notes, numpy overflow) are not printed
mkdir -p logs
echo "run_all.sh started $(date)"
echo "python $(python -c 'import sys; print(sys.version.split()[0])'), python-flint $(python -c 'import flint; print(flint.__version__)'), numpy $(python -c 'import numpy; print(numpy.__version__)'), cvxpy $(python -c 'import cvxpy; print(cvxpy.__version__)' 2>/dev/null), sympy $(python -c 'import sympy; print(sympy.__version__)')"
allok=1
step () {
  name="$1"; lim="$2"; want="$3"; shift 3
  t0=$(date +%s)
  bash ./guard.sh "logs/run_${name}.log" "$lim" "${name}.py" "$@"
  st=$?
  t1=$(date +%s)
  last=$(tail -n 1 "logs/run_${name}.log" | tr -d '\r')
  if [ "$st" = 0 ] && [ "$last" = "$want" ]; then res=PASS; else res=FAIL; allok=0; fi
  echo "step ${name}: ${res} (exit status ${st}, $((t1 - t0)) s, last line: ${last})"
}
step check_j1_ylib 900 PASS                # the Lemma B builder j1lib.py equals ylib.Builder on six |T| = 1 blocks
step check_q 1800 "ALL PASS"               # Proposition Q and Lemma K, exact; symbolic identity for |U| = 2, 3
step check_rp 900 "ALL PASS"               # Proposition RP, exact
step cert_pairs 1800 "ALL PASS" 1          # exact pair certificates (Section 4) and negative controls
step num_pairs 2400 PASS 40                # floating-point survey for Conjecture TC (numerical, no claim rests on it)
step explore_structure 1800 "numerical exploration done"   # floating-point remarks of Section 5 (numerical)
step explore_rq 1800 "numerical exploration done"          # 50-digit remark on the (R, Q) sub-block (numerical)
echo "run_all.sh finished $(date)"
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
