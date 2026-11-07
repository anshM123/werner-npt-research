#!/usr/bin/env bash
# Fresh rerun of the steps of run_st4.sh from a clean copy (no slices, no results of stage 4), with the stage-3
# results res/{anti,sym1}_{MT,MU}.pkl, res/det_MT.pkl, res/det1_MT.pkl taken from the fresh rerun of stage 3.  Only the
# families the claims rest on are run: step 4 on the two proved D3 families, step 5 on the 11 proved MU families; the
# two deferred D3 families keep their STOP files.  Usage: NW=2 bash run_verify_st4.sh > logs/verify_st4.out
cd "$(dirname "$0")"
NW=${NW:-2}
echo "run_verify_st4.sh started $(date)"
NW=$NW bash run_det1_large.sh > logs/run_det1_large.out; cat logs/run_det1_large.out
NW=$NW bash run_sd.sh lists/MU_sym_2d_fams.txt > logs/run_sd_MU2d.out; echo "sd 2d: $(grep -c ' / PASS$' logs/run_sd_MU2d.out) PASS of $(grep -c . lists/MU_sym_2d_fams.txt)"
NW=$NW bash run_d2.sh lists/MT_sym_3d.txt > logs/run_d2_MT.out; echo "d2: $(grep -c '^final .*: PASS$' logs/run_d2_MT.out) PASS of 10"
NW=$NW bash run_d3.sh lists/MT_sym_3e_proved.txt > logs/run_d3_MT.out; echo "d3: $(grep -c '^final .*: PASS$' logs/run_d3_MT.out) PASS of 2"
NW=$NW bash run_sd.sh lists/MU_sym_2e_cheap_proved.txt > logs/run_sd_MU2e.out; echo "sd 2e: $(grep -c ' / PASS$' logs/run_sd_MU2e.out) PASS of 11"
bash guard.sh logs/run_q3cert.log 5400 q3cert.py res3/q3.pkl @lists/q3_list.txt; echo "q3cert: $(tail -n 1 logs/run_q3cert.log)"
bash guard.sh logs/run_negdet.log 3600 negdet.py $(tr -d '\r' < lists/neg_list.txt); echo "negdet: $(tail -n 1 logs/run_negdet.log) (expected FAIL)"
bash guard.sh logs/run_coverage_st4.log 3600 coverage_st3.py 120 29 30 res/*.pkl res3/q3.pkl; echo "coverage: $(tail -n 1 logs/run_coverage_st4.log)"
echo "run_verify_st4.sh done $(date)"
