#!/usr/bin/env bash
# Second check of Corollaries TV and M1: every torus group of (R_inf)_{k-1} at two weight values (N copies of amplitude
# A), from Lemma B with exact positive-definiteness certificates.  One Python process at a time under guard.sh.
cd "$(dirname "$0")"
mkdir -p logs
for a in "10 3 3/2 9" "10 5 2/5 9" "11 4 3 10" "11 1 1/4 10" "12 5 2/3 11" "12 1 7/3 11" "13 6 1/2 12"; do
  set -- $a
  name="markall2_k$1_n$2_$(echo $3 | tr '/' 'o')"
  bash guard.sh logs/$name.log 14400 markall2.py $1 $2 $3 $4
  echo "$name: $(tail -n 1 logs/$name.log) $(grep '^total' logs/$name.log)"
done
