#!/bin/sh
# Theorem U2 (NOTES.md): G > 0 for the type I blocks of (R_inf)_{k-1} with |U| = 2, every m, all weights.
# Reruns every step, one Python process at a time (each well under 1 GB), prints PASS/FAIL per script and ALL PASS.
# The closed forms (closed_u2.pkl) and the polynomial data (build_A.pkl, build_B.pkl) are regenerated here.
# Needs Python 3 with numpy and python-flint (0.9 used here).  About 1 minute.
cd "$(dirname "$0")"
mkdir -p logs
: > logs/run_all.log
fail=0
run() {
  name=$1; shift
  start=$(date +%s)
  python -u ${name}.py "$@" > logs/${name}.log 2>&1
  last=$(tail -n 1 logs/${name}.log)
  end=$(date +%s)
  echo "$name $*: $last ($((end - start)) s)" | tee -a logs/run_all.log
  [ "$last" = "PASS" ] || fail=1
}
run test_bern            # Bernstein module unit tests (copy of the repo module)
run check_cf             # closed form cf.py = Lemma B builder (exact) and = direct construction (float)
run closed_u2            # exact symbolic Schur complement (u = 2, m = 3) -> closed_u2.pkl; evenness; symmetry
run check_closed         # closed forms = direct exact elimination for m = 1..5 (115 blocks); norm formulas (float)
run build_B              # (B_2): minors as polynomials -> build_B.pkl
run prove_B 2/5          # (B_2): denominator signs, blow-up charts, Bernstein certificates, one-point check
run build_A 2/5          # (A_2): the polynomial FA -> build_A.pkl
run prove_A 2/5          # (A_2): denominator signs, triangle map, Bernstein certificate
run controls             # negative controls H = 1/4 (A) and H = 4/5 (B) with exact witnesses
if [ $fail -eq 0 ]; then echo "ALL PASS" | tee -a logs/run_all.log; else echo "SOME FAIL" | tee -a logs/run_all.log; fi
