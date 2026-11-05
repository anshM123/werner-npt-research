#!/bin/sh
# G > 0 for type I blocks of (R_inf)_{k-1}: the case |U| = 3, |T| = 1 (k = 4), plus the same certificate for small (u, m).
# Reruns every step, one Python process at a time, prints PASS/FAIL per script and ALL PASS / SOME FAIL.
# All closed forms and polynomial data are regenerated here.  Needs Python 3 with numpy and python-flint (0.9 used).
# FULL=1 also attempts (u, m) = (4, 1), which did not finish in our time budget (see NOTES.md, Section 5).
cd "$(dirname "$0")"
mkdir -p logs
: > logs/run_all.log
fail=0
run() {
  name=$1; shift
  tag=$(echo "$name $*" | sed 's/ *$//' | tr ' ' '_')
  start=$(date +%s)
  python -u ${name}.py "$@" > logs/${tag}.log 2>&1
  last=$(tail -n 1 logs/${tag}.log | tr -d '\r')
  end=$(date +%s)
  echo "$name $*: $last ($((end - start)) s)" | tee -a logs/run_all.log
  [ "$last" = "PASS" ] || fail=1
}
run test_bern                 # unit tests of the Bernstein module (used only as a fallback in prove_m1.py)
run check_cf                  # closed form cf.py = Lemma B builder (exact) and = direct construction (float)
# ---- the theorem: u = 3, m = 1 (k = 4)
run closed_m1 3               # exact Schur complement G (4 x 4) over Q(s0..s3) -> closed_m1_u3.pkl  (about 1 min)
run check_closed_m1 3         # closed form = direct exact elimination at 12 rational weight vectors, two gauges
run build_m1                  # leading minors of den*G_z with the denominator factors stripped -> build_m1.pkl
run prove_m1                  # signs: den factors and all leading-minor cores have positive coefficients => G > 0
run controls_m1               # one-point LDL check; negative controls with exact witnesses
# ---- the same certificate for small (u, m): an independent re-proof of known cases, and the pattern
for um in "1 1" "2 1" "1 2" "2 2" "1 3" "3 1"; do
  run closed_gen $um          # exact Schur complement for (u, m) -> closed_u{u}_m{m}.pkl
  run check_closed_gen $um    # closed form = direct exact elimination at 6 rational points, two gauges
  run minors_pos $um          # leading minors, factor stripping, coefficient signs => G > 0 for (u, m)
done
run check_zpattern            # exact: for m = 1, u <= 3, G_z is a symmetric Z-matrix (Section 4b)
run check_mvec                # exact: for m = 1, u <= 3, G_z v > 0 for an explicit v (second proof via an M-matrix, Section 4b)
if [ "${FULL:-0}" = "1" ]; then
  run closed_gen 4 1
  run check_closed_gen 4 1 3
  run minors_pos 4 1
fi
# ---- numerical evidence (floating point; no claim rests on it)
run sweep_u                   # Gershgorin constants, lambda_min, perp block, merged comparison, u <= 7, m <= 4
run ex_mnum                   # m = 1, u <= 7: Z-matrix pattern and the window for v = (1, beta, ..., beta) (floating point)
run ex_mvec                   # m = 1, u <= 8: the explicit v of check_mvec.py still gives G_z v > 0 (floating point)
if [ $fail -eq 0 ]; then echo "ALL PASS" | tee -a logs/run_all.log; else echo "SOME FAIL" | tee -a logs/run_all.log; fi
