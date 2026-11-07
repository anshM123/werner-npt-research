#!/usr/bin/env bash
# development helper (not part of any proof): hands the process slots freed by the D3 run and by q3cert to more
# workers of the MU Proposition D run, keeping at most six Python processes of these runs at once
cd "$(dirname "$0")/.."
until ls logs/d3_final_*.log > /dev/null 2>&1; do sleep 5; done
for nm in $(sed 's/held://' logs/dev/held_claims3.txt); do rmdir "claims/sd_$nm" 2>/dev/null; done
echo "$(date) d3 workers done: MU worker A"
JOIN=1 NW=1 bash run_sd.sh lists/MU_sym_2d_fams.txt > logs/run_sd_MU2d_part4.out 2>&1 &
until grep -q "run_d3.sh .* done" logs/run_d3_MT.out; do sleep 5; done
echo "$(date) d3 finals done: MU worker B"
JOIN=1 NW=1 bash run_sd.sh lists/MU_sym_2d_fams.txt > logs/run_sd_MU2d_part5.out 2>&1 &
until [ -n "$(tail -n 1 logs/run_q3cert.log | grep -E '^(PASS|FAIL)')" ] || grep -q "guard:" logs/run_q3cert.log; do sleep 5; done
echo "$(date) q3cert done: MU worker C"
JOIN=1 NW=1 bash run_sd.sh lists/MU_sym_2d_fams.txt > logs/run_sd_MU2d_part6.out 2>&1 &
wait
echo "$(date) all MU workers done"
