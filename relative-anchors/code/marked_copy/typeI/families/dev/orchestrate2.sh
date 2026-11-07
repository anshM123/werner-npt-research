#!/usr/bin/env bash
# development helper (not part of any proof): after the MU Proposition D workers, the cheap MU families of
# Proposition D1 (three workers); after the D2 run, the negative controls and two more workers; at most six Python
# processes of these runs at once
cd "$(dirname "$0")/.."
until grep -q "all MU workers done" logs/dev/orchestrate1.out; do sleep 5; done
echo "$(date) MU Proposition D done: three workers for lists/MU_sym_2e_cheap.txt"
NW=3 bash run_sd.sh lists/MU_sym_2e_cheap.txt > logs/run_sd_MU2e_cheap.out 2>&1 &
until grep -q "run_d2.sh .* done" logs/run_d2_MT.out; do sleep 5; done
echo "$(date) D2 done: negative controls and two more workers"
JOIN=1 NW=2 bash run_sd.sh lists/MU_sym_2e_cheap.txt > logs/run_sd_MU2e_cheap_join.out 2>&1 &
bash guard.sh logs/run_negdet.log 3600 negdet.py $(tr -d '\r' < lists/neg_list.txt)
echo "$(date) negdet: $(tail -n 1 logs/run_negdet.log)"
wait
echo "$(date) all done"
