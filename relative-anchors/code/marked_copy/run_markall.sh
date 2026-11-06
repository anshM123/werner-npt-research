#!/usr/bin/env bash
# Second check of Theorem I1 / Corollary I2: every torus group of (R_inf)_j, 3 <= j <= k - 1, at weights all equal
# but one, from Lemma B with exact positive-definiteness certificates.  One Python process at a time under guard.sh.
cd "$(dirname "$0")"
mkdir -p logs
for a in "8 2/3 3,4,5,6,7" "8 3 3,4,5,6,7" "9 1/2 3,4,5,6,7,8" "9 5/2 3,4,5,6,7,8" "10 3/2 3,8,9"; do
  set -- $a
  name="markall_k$1_$(echo $2 | tr '/' 'o')"
  bash guard.sh logs/$name.log 14400 markall.py $1 $2 $3
  echo "$name: $(tail -n 1 logs/$name.log) $(grep '^total' logs/$name.log)"
done
