#!/usr/bin/env bash
# Reruns every step of NOTES.md, one Python process at a time, each under guard.sh (memory limit 2 500 000 KB and a
# time limit).  Prints PASS or FAIL per step and ALL PASS at the end.
# Usage: bash run_all.sh > logs/run_all.out 2>&1      (one log logs/run_<step>.log per step)
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
# 1. closed forms (Sections 2-3): u1 = 1 (general T-weights; T-weights in {s, 1}), and u1, u2 >= 2 (T-weights in {s, 1})
step tl_inv_one 1500 PASS tl_m.py inv one
for md in T11 T12 T22; do
  step tl_inv_one_$md 900 PASS tl_m.py inv one $md
  step tl_std2_one_$md 900 PASS tl_m.py std2 one $md
  step tl_std2_onetwo_$md 900 PASS tl_m.py std2 onetwo $md
  step tl_inv_gen_$md 1500 PASS tl_m.py inv gen $md
  step tl_std1_gen_$md 900 PASS tl_m.py std1 gen $md
  step tl_std1_two_$md 900 PASS tl_m.py std1 two $md
done
# 2. exact validation of every closed form against direct elimination of whole blocks (Section 3)
step check_tl_one 1500 PASS check_tl.py tl_inv_one.pkl tl_inv_one_T11.pkl tl_inv_one_T12.pkl tl_inv_one_T22.pkl tl_std2_one_T11.pkl tl_std2_one_T12.pkl tl_std2_one_T22.pkl tl_std2_onetwo_T11.pkl tl_std2_onetwo_T12.pkl tl_std2_onetwo_T22.pkl
step check_tl_gen 2400 PASS check_tl.py tl_inv_gen_T11.pkl tl_inv_gen_T12.pkl tl_inv_gen_T22.pkl tl_std1_gen_T11.pkl tl_std1_gen_T12.pkl tl_std1_gen_T22.pkl tl_std1_two_T11.pkl tl_std1_two_T12.pkl tl_std1_two_T22.pkl
# 3. Theorem M1: one marked copy in U, all other weights equal (exact S_m splitting, Section 4)
step prove_eqT 900 PASS prove_eqT.py
# 4. Theorem TV: two-valued weight vectors (weighted block Gershgorin, H = 1/2, Section 5)
step prove_gersh_one 1800 PASS prove_gersh.py one
step prove_gersh_one_H2_5 1800 PASS prove_gersh.py one 2 5
# u1, u2 >= 2, both symbolic: one process per group of conditions (memory).  (B_inv) only its first two pivots
# (the third one exceeds the memory limit with both sizes symbolic; NOTES.md, Section 5)
for tag in A_ss_rho A_ss_zz A_s1_rho A_s1_zz A_11_rho A_11_zz B_s12 B_112 std; do
  step prove_gersh_gen_$tag 1800 PASS prove_gersh.py gen 1 2 only=$tag
done
# (B_inv) completely (three pivots) for u1 = N fixed, u2 >= 2 symbolic, N = 2..NMAX; with the exchange of the
# levels this gives (B_inv) whenever min(u1, u2) <= NMAX
NMAX=${NMAX:-16}
for n in $(seq 2 $NMAX); do
  step tl_inv_f${n}_T11 900 PASS tl_m.py inv f$n T11
  step tl_inv_f${n}_T22 900 PASS tl_m.py inv f$n T22
  step prove_gersh_f$n 900 PASS prove_gersh.py f$n
done
# the fixed-u1 forms: against direct elimination (u1 = 2, 3, 5) and against the symbolic forms at u1 = N (all N)
step check_tl_fN 2400 PASS check_tl.py tl_inv_f2_T11.pkl tl_inv_f2_T22.pkl tl_inv_f3_T11.pkl tl_inv_f3_T22.pkl tl_inv_f5_T11.pkl tl_inv_f5_T22.pkl
step check_fN 900 PASS check_fN.py $NMAX
# 5. end-to-end checks and negative controls (Section 6)
step check_e2e2 2400 PASS check_e2e2.py
step controls2 2400 PASS controls2.py
# 6. floating-point exploration (no claim rests on it)
step ex_numframe 900 "exploration done" ex_numframe.py
step ex_twoval 900 "exploration done" ex_twoval.py
step ex_factor_one 900 "exploration done" ex_factor_one.py
step ex_signs_one 900 "exploration done" ex_signs_one.py
echo "run_all.sh finished $(date)"
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
