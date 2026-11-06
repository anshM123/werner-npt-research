#!/usr/bin/env bash
# batch 1: every anti family, then the sym families with one symbolic size; one Python process at a time under guard.sh
cd "$(dirname "$0")"
for lst in MT_anti_1 MT_anti_2 MT_anti_3 MT_anti_4 MU_anti_1 MU_anti_2 MU_anti_4 MU_anti_3 MT_sym_1 MU_sym_1; do
  t0=$(date +%s)
  bash guard.sh logs/fam_$lst.log 3600 tprove.py res/fam_$lst.pkl @lists/$lst.txt --tlimit=600
  st=$?
  echo "$lst exit $st $(tail -n 2 logs/fam_$lst.log | head -n 1 | tr -d '\r') $(tail -n 1 logs/fam_$lst.log | tr -d '\r') $(( $(date +%s) - t0 ))s"
done
echo BATCH1 DONE
