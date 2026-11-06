#!/usr/bin/env bash
# A check of Theorems AT, S1, SD, S2 written separately from their certificates (st3check.py): at three sampled members
# of every certified family with 30 <= k <= 90 (two seeds), the numeric reduced form of tlib.py at a rational sigma is
# positive semidefinite with kernel exactly the span of its known kernel vectors.  Run after run_all.sh and
# run_extra.sh (it reads res/*.pkl).  One process at a time under guard.sh.
cd "$(dirname "$0")"
mkdir -p logs
for seed in 11 12; do
  bash guard.sh logs/st3check_$seed.log 14400 st3check.py 3 $seed 90 res/anti_MT.pkl res/anti_MU.pkl res/sym1_MT.pkl res/sym1_MU.pkl res/det_MT.pkl res/det1_MT.pkl
  echo "seed $seed: $(tail -n 2 logs/st3check_$seed.log | tr '\n' ' ')"
done
