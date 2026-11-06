#!/usr/bin/env bash
# quick tests of the code paths of run_extra.sh (one Python process at a time, under guard.sh); not part of any proof
cd "$(dirname "$0")"
mkdir -p logs/smoke
bash guard.sh logs/smoke/s1.log 300 structure.py --dens MT/anti/1,1,n,n MU/anti/2,2,n,n; tail -2 logs/smoke/s1.log
bash guard.sh logs/smoke/s2.log 300 tdet.py prove logs/smoke/det.pkl MT 1,1; tail -4 logs/smoke/s2.log
bash guard.sh logs/smoke/s3.log 300 tdet.py measure MT 1 3; tail -4 logs/smoke/s3.log
bash guard.sh logs/smoke/s4.log 300 tdet.py measure_plain MT/sym/n,1,1,0; tail -3 logs/smoke/s4.log
bash guard.sh logs/smoke/s5.log 300 probe.py MT/sym/n,1,1,0 20; tail -2 logs/smoke/s5.log
bash guard.sh logs/smoke/s6.log 300 summary.py res/anti_MT.pkl logs/smoke/det.pkl; tail -3 logs/smoke/s6.log
echo SMOKE DONE
