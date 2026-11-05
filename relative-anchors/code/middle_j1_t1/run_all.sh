#!/usr/bin/env bash
# Reruns every check of NOTES.md, one Python process at a time, each under guard.sh (killed above 2.5 GB or after its
# time limit).  Prints PASS/FAIL per step and ALL PASS at the end.  Per-step logs: logs/run_<step>.log; the summary is
# saved as logs/run_all.out (bash run_all.sh > logs/run_all.out).  About 12 minutes (6 of them in cert_F_l2_22_inv).
# LONG=1 also starts the attempts that did not complete in our runs: the k = 5 all-weights certificate and the
# three-level family (NOTES.md, Sections 3b and 4).
# The floating-point explorations (explore_*.py) carry no claim and are not run here; their logs are logs/explore_*.log.
# Needs Python 3 with python-flint (0.9 used); numpy and mpmath only for the explorations.
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
# the builders: j1lib.py (Lemma B) against ylib.Builder; Lemma T1.3 against j1lib (copied check, unchanged)
step check_j1_ylib 600 check_j1_ylib.py
step check_t1_closed 900 check_t1_closed.py
# Lemma P and the closed form of Phi, exact against j1lib
step check_phi_closed 900 check_phi_closed.py
# the level reductions of lev.py against j1lib (one, two and three levels; B-part scaled by 1 and 1 + 1/1000)
step check_lev 900 check_lev.py
# Theorem EQ: the separate first implementation, its certificate, and the comparison with lev.py
step equ_sym 300 equ_sym.py
step cert_equ 300 cert_equ.py
step cmp_impl 300 cmp_impl.py
# Theorem EQ: certificate 1 (Phi) and certificate 2 (gauge-fixed block, P_t deleted)
step cert_F_l1 600 cert_lev.py F 1 2
step cert_G_l1 600 cert_lev.py G 1 2
# negative controls: the certificate must fail with the B-part scaled by 1 + 1/1000; explicit witnesses
./guard.sh logs/run_cert_G_l1_neg.log 600 -u cert_lev.py G 1 2 1001/1000
if tail -1 logs/run_cert_G_l1_neg.log | tr -d '\r' | grep -q '^FAIL$'; then
  echo "PASS  cert_G_l1_neg (the certificate fails for the scaled B-part, as it must)"
else
  echo "FAIL  cert_G_l1_neg (the certificate did not fail for the scaled B-part)"; fail=1
fi
step neg_lev 900 neg_lev.py
# Theorem 2L: u1, u2 >= 2 (invariant part, about 6 minutes); standard parts; u1 = 1 exactly; u1 = u2 = 1 (k = 3)
step cert_F_l2_22_inv 1800 cert_lev.py F 2 2 2 only=inv
step cert_F_l2_21_std 600 cert_lev.py F 2 2 1 only=std
step cert_F_l2_12_std 600 cert_lev.py F 2 1 2 only=std
step cert_F_l2_12_fixed_inv 900 cert_lev.py F 2 1 2 fixed=1,0 only=inv
step cert_F_l2_11_fixed 600 cert_lev.py F 2 1 1 fixed=1,1 only=inv
# fixed k, all weights, through Phi' (Section 4): k = 3, 4
step cert_pivots_k3 300 cert_pivots.py 3
step cert_pivots_k4 1200 cert_pivots.py 4
step bareiss_phic_k4 900 bigstack.py bareiss_phic.py 4
if [ "${LONG:-0}" = "1" ]; then
  # not completed in our runs (NOTES.md, Sections 3b and 4): memory above 2.5 GB, or a segmentation fault in the
  # polynomial library
  step bareiss_phic_k5 5400 bigstack.py bareiss_phic.py 5
  step bareiss_l3_inv 3600 bigstack.py bareiss_phic.py levels 3 1,1,2 1,1,0 inv F
  step pivots_l3_inv 3600 bigstack.py pivots_lev.py 3 1,1,2 1,1,0 inv
fi
if [ $fail -eq 0 ]; then echo "ALL PASS"; else echo "SOME STEPS FAILED"; fi
