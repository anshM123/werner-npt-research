#!/usr/bin/env bash
# A check of Proposition B3 by the other method: the fixed-u1 certificate of Section 5 (prove_gersh.py fN: all three
# pivots of (B_inv), u2 >= 2 symbolic, by elimination with coefficient signs) at values of N beyond 16, inside and
# outside the interpolation grid u1 = 2..78.  One Python process at a time under guard.sh.
cd "$(dirname "$0")"
mkdir -p logs
for n in 17 18 20 24 30 40 50 77 100 200 500; do
  for md in T11 T22; do
    bash guard.sh logs/tl_inv_f${n}_$md.log 1800 tl_m.py inv f$n $md
    echo "tl_inv_f${n}_$md: $(tail -n 1 logs/tl_inv_f${n}_$md.log | tr -d '\r')"
  done
  bash guard.sh logs/prove_gersh_f$n.log 1800 prove_gersh.py f$n
  echo "prove_gersh_f$n: $(tail -n 1 logs/prove_gersh_f$n.log | tr -d '\r')"
done
