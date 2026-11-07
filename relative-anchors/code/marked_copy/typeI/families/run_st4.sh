#!/usr/bin/env bash
# The steps of NOTES.md, Sections 8 to 10 (after run_all.sh and run_extra.sh of Sections 1 to 7, whose results
# res/*.pkl are those of stage 3).  Each step is a run script with parallel workers, each worker one Python process
# under guard.sh (1.2 GB); keep NW + (other jobs) within the machine's limits.  Every step is resumable: finished
# families (log ending with PASS, or the marker final_ok in their folder) and saved slices are skipped.
# Usage: NW=4 bash run_st4.sh > logs/run_st4.out
# (In the session of 6 November the same steps were run with up to six processes at once, in the order 1, 3, 4, 2,
# 5, 6, 7, 8 below, partly overlapping; the logs named here are those of that session.)
cd "$(dirname "$0")"
NW=${NW:-4}
mkdir -p logs res res3
t_start=$(date +%s)
echo "run_st4.sh started $(date)"
# 1. the 19 MT sym families with two symbolic sizes, the fixed one of i, tau' equal to 1 (Proposition D1, Section 8.2)
[ -f res/det1_MT_large_w0.pkl ] || NW=3 bash run_det1_large.sh > logs/run_det1_large.out
# 2. the 24 MU sym families [q0 = a, t = c, i >= 2, tau' >= 2] (Proposition D with the split, Section 8.3)
NW=$NW bash run_sd.sh lists/MU_sym_2d_fams.txt > logs/run_sd_MU2d.out
# 3. the 10 MT sym families with i, tau' and one of q0, t symbolic (Proposition D2, Section 8.4)
NW=$NW bash run_d2.sh lists/MT_sym_3d.txt > logs/run_d2_MT.out
# 4. the other 4 MT sym families with three symbolic sizes (Proposition D3, Section 8.5); the folders
#    d3/MT_sym_n_n_n_1 and d3/MT_sym_n_n_1_n hold a file STOP (deferred): remove it to run them
NW=$NW bash run_d3.sh lists/MT_sym_3e.txt > logs/run_d3_MT.out
# 5. the 14 cheap MU sym families with two symbolic sizes, not i and tau' (Proposition D1 with the split, Section 8.5b)
NW=$NW bash run_sd.sh lists/MU_sym_2e_cheap.txt > logs/run_sd_MU2e_cheap.out
# 6. the three-row parts of every family that may be open (Section 8.6)
bash guard.sh logs/run_q3cert.log 5400 q3cert.py res3/q3.pkl @lists/q3_list.txt
echo "q3cert: $(tail -n 1 logs/run_q3cert.log)"
# 7. negative controls (Section 10)
bash guard.sh logs/run_negdet.log 3600 negdet.py $(tr -d '\r' < lists/neg_list.txt)
echo "negdet: $(tail -n 1 logs/run_negdet.log) (expected FAIL)"
# 8. coverage with every result (Section 9)
bash guard.sh logs/run_coverage_st4.log 3600 coverage_st3.py 120 29 30 res/*.pkl res3/q3.pkl
echo "coverage: $(tail -n 1 logs/run_coverage_st4.log)"
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
