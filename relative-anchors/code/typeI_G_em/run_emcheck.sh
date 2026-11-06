#!/usr/bin/env bash
# Second check of Corollary EM: whole (R_inf)_(k-1) from Lemma B at weights equal on |U| >= 3 copies, k = 5..8.
cd "$(dirname "$0")"
mkdir -p logs
for a in "5 10 51" "6 8 61" "7 4 71" "8 2 81"; do
  set -- $a
  bash guard.sh logs/emcheck_k$1.log 7200 emcheck.py $1 $2 $3
  echo "k=$1: $(tail -n 1 logs/emcheck_k$1.log) $(grep '^total' logs/emcheck_k$1.log)"
done
