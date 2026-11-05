#!/usr/bin/env bash
# Reruns every check of NOTES.md, one Python process at a time, each under guard.sh (killed above 2.5 GB or after its
# time limit).  Prints PASS/FAIL per step and ALL PASS at the end.  Per-step logs: logs/run_<step>.log; the summary is
# saved as logs/run_all.out (bash run_all.sh > logs/run_all.out).  About three minutes.
# Needs Python 3 with python-flint (0.9 used).
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
# a control step: the certificate with the B-part scaled by 1 + 1/1000 must FAIL
negstep() {
  name="$1"; lim="$2"; shift 2
  ./guard.sh "logs/run_${name}.log" "$lim" -u "$@"
  st=$?
  if [ $st -eq 0 ] && tail -1 "logs/run_${name}.log" | tr -d '\r' | grep -q '^FAIL$'; then
    echo "PASS  $name (the certificate fails for the scaled B-part, as it must)"
  else
    echo "FAIL  $name (exit $st; the certificate did not fail for the scaled B-part; see logs/run_${name}.log)"
    fail=1
  fi
}
# 1. the builders: j1lib (Lemma B) against ylib (definitions), general weights, k = 4, 5, 6 (copied check, unchanged)
step check_ylib 900 check_ylib.py
# 2. the closed forms of l1sym.Model entry by entry against j1lib (B-part scaled by 1 and 1 + 1/1000)
step check_entries 900 check_entries.py
# 3. the reduced forms (pattern sums) against explicit projections of the j1lib block, 14 size pairs, both scalings
step check_reduce 900 check_reduce.py
# 3b. the same at two larger size pairs, (5, 4) and (6, 3) (k = 9), irreps on both levels
step check_reduce_large 1800 check_reduce_large.py
# 4. the reduced forms over Q(q) at fixed sizes against ylib over Q(q), identically in q
step check_ylib_sym 900 check_ylib_sym.py
# 5. the symbolic forms of the certificates at rational points against j1lib (both orders of the levels), and at
#    larger sizes against the reduction in rational arithmetic
step check_sym_points 1800 check_sym_points.py
# 6. bookkeeping: the certified regimes cover every part at every (u0, u1) with k >= 4; dimension counts
step coverage 300 coverage.py
# 7. the certificates (Theorem L1-2): every regime
for r in "inv 3+ 1" "inv 2 2" "inv 3+ 2" "inv 3+ 3+" "inv 4+ 0" \
         "std 3+ 1" "std 2 2" "std 2 3+" "std 3+ 2" "std 3+ 3+" "std 4+ 0" \
         "two 4+ 1" "two 4+ 2" "two 4+ 3+" "two 4+ 0" \
         "anti 3+ 1" "anti 3+ 2" "anti 3+ 3+" "anti 4+ 0" \
         "mix 2 2" "mix 2 3+" "mix 3+ 3+"; do
  n=$(echo "$r" | tr ' +' '_p')
  step "cert_$n" 1800 cert.py $r
done
# 8. negative controls: the same certificates with the B-part scaled by 1 + 1/1000 must fail (the parts that carry
#    the kernel, and 'two' at u1 = 1); explicit rational points where the scaled gauge-fixed block is not PD
for r in "inv 3+ 1" "inv 2 2" "inv 3+ 2" "inv 3+ 3+" "std 3+ 1" "std 2 2" "std 2 3+" "std 3+ 2" "std 3+ 3+" \
         "two 4+ 1"; do
  n=$(echo "$r" | tr ' +' '_p')
  negstep "neg_cert_$n" 1800 cert.py $r 1001/1000
done
step neg 900 neg.py
# observation (row 9 of the status table, not used in any proof): these certificates still pass with the scaled B-part
for r in "two 4+ 2" "two 4+ 3+" "anti 3+ 1" "anti 3+ 2" "anti 3+ 3+" "mix 2 2" "mix 2 3+" "mix 3+ 3+"; do
  n=$(echo "$r" | tr ' +' '_p')
  step "obs_scaled_$n" 1800 cert.py $r 1001/1000
done
# 9. the conclusion on the whole block at 54 two-level rational weight vectors, k = 4..7 (exact checks)
step check_block 1800 check_block.py
if [ $fail -eq 0 ]; then echo "ALL PASS"; else echo "SOME STEPS FAILED"; fi
