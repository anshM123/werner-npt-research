#!/usr/bin/env bash
# development helper (not part of any proof): when the D2 run is done, the negative controls and two workers for the
# cheap MU families of Proposition D1; after the controls one more worker; when the MU Proposition D run is done three
# more; at the end the coverage.  At most six Python processes of these runs at once (the MU Proposition D run has at
# most three workers).
cd "$(dirname "$0")/.."
until grep -q "run_d2.sh .* done" logs/run_d2_MT.out; do sleep 5; done
echo "$(date) D2 done: negdet and two workers for the cheap MU families"
NW=2 bash run_sd.sh lists/MU_sym_2e_cheap.txt > logs/run_sd_MU2e_cheap.out 2>&1 &
bash guard.sh logs/run_negdet.log 3600 negdet.py $(tr -d '\r' < lists/neg_list.txt)
echo "$(date) negdet: $(tail -n 1 logs/run_negdet.log); one more worker"
JOIN=1 NW=1 bash run_sd.sh lists/MU_sym_2e_cheap.txt > logs/run_sd_MU2e_cheap_j1.out 2>&1 &
until grep -q "all MU workers done" logs/dev/orchestrate1.out; do sleep 5; done
echo "$(date) MU Proposition D done: three more workers"
JOIN=1 NW=3 bash run_sd.sh lists/MU_sym_2e_cheap.txt > logs/run_sd_MU2e_cheap_j2.out 2>&1 &
wait
echo "$(date) cheap MU families done: coverage"
bash guard.sh logs/run_coverage_st4.log 3600 coverage_st3.py 120 29 30 res/*.pkl res3/q3.pkl
echo "$(date) coverage: $(tail -n 1 logs/run_coverage_st4.log)"
echo "$(date) all done"
