#!/usr/bin/env bash
# The check written separately (st3check.py) on the families certified in stage 4 (Theorems S3, SD2, MD, the two
# families of Proposition D3 and the 11 of Section 8.5b): three sampled members of every family with 30 <= k <= 90 for
# two seeds, and ten members with k <= 150 for a third.  Run after run_verify_st4.sh (it reads res/*.pkl).
cd "$(dirname "$0")"
mkdir -p logs
FILES=$(ls res/det1_MT_large_w*.pkl res/sd_*.pkl res/det2_*.pkl res/det3_*.pkl)
for seed in 21 22; do
  bash guard.sh logs/st4check_$seed.log 7200 st3check.py 3 $seed 90 $FILES
  echo "seed $seed: $(tail -n 2 logs/st4check_$seed.log | tr '\n' ' ')"
done
bash guard.sh logs/st4check_wide.log 7200 st3check.py 10 23 150 $FILES
echo "wide: $(tail -n 2 logs/st4check_wide.log | tr '\n' ' ')"
