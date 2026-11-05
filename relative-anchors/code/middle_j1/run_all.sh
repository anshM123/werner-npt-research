#!/usr/bin/env bash
# Reruns every check of NOTES.md, one Python process at a time, each under guard.sh (killed above 2.8 GB or after its
# time limit).  Prints PASS/FAIL per step and ALL PASS at the end.  Logs: logs/run_<step>.log; the summary is saved
# as logs/run_all.out (bash run_all.sh | tee logs/run_all.out).  About 55 minutes, 52 of them in cert_t3j2.py.
# LONG=1 also starts the k = 5 leading-minor certificate for |T| = 1 (bareiss_t1.py 5, through bigstack.py) with a
# one-hour limit; in our runs it did not finish (NOTES.md, Section 3), so it is off by default.
# The floating-point explorations (explore_*.py) carry no claim and are not run here; their logs are logs/explore_*.log.
# Needs Python 3 with numpy, sympy and python-flint (0.9 used).
cd "$(dirname "$0")"
mkdir -p logs
fail=0
step() {
  name="$1"; lim="$2"; shift 2
  start=$(date +%s)
  ./guard.sh "logs/run_${name}.log" "$lim" -u "$@"
  st=$?
  end=$(date +%s)
  if [ $st -eq 0 ] && tail -1 "logs/run_${name}.log" | tr -d '\r' | grep -q '^PASS$'; then
    echo "PASS  $name ($((end - start)) s)"
  else
    echo "FAIL  $name (exit $st; see logs/run_${name}.log)"
    fail=1
  fi
}
# Lemma B builder against ylib (k = 4, 5, 6), block sizes and kernels
step check_ylib 900 check_ylib.py
# Theorem D1 (|T| = 2, every k): certificate of the eight reduced inequalities, and the exact checks
step cert_t2 600 cert_t2.py
step check_t2 900 check_t2.py
# |T| = 1: closed form and inverse-free reduction (exact); Conjecture L exact checks k = 3..10
step check_t1_closed 900 check_t1_closed.py
step check_t1 2400 check_t1.py 10 15
# |T| = 1: Conjecture L proved at k = 3, 4 (symbolic), and the leading-minor certificate at k = 4
step sym_t1_k3 300 sym_t1.py 3
step sym_t1_k4 600 sym_t1.py 4
step bareiss_t1_k4 600 bareiss_t1.py 4
if [ "${LONG:-0}" = "1" ]; then
  step bareiss_t1_k5 3600 bigstack.py bareiss_t1.py 5
fi
# letter-1 block: vacuum reduction and inverse-free form (exact); Conjecture L1 exact checks k = 4, 5, 6; control
step check_l1_closed 900 check_l1_closed.py
step check_l1 1800 check_l1.py 6 4
# Theorem D2 ((R_inf)_2, |T| = 3, every k): the 8 + 64 vertex certificates (about 52 minutes) and the exact checks
step cert_t3j2 7200 cert_t3j2.py   # took 3103 s in the recorded run; the limit leaves room
step check_t3j2 1200 check_t3j2.py
# the formulas certified by cert_t2.py and cert_t3j2.py, evaluated exactly against the Lemma B blocks
step check_cert_formulas 900 check_cert_formulas.py
if [ $fail -eq 0 ]; then echo "ALL PASS"; else echo "SOME STEPS FAILED"; fi
