#!/usr/bin/env bash
# Stage 2 of NOTES Section 6c: the letter-1 block of (R_inf)_2 at two weight values with BOTH level sizes symbolic,
# u0 = 14 + x0, u1 = 14 + x1 (x0, x1 >= 0 real): the certificate of cert.py for every family, in (q, x0, x1).
# One Python process at a time under guard.sh (2.5 GB, 10 hours per family).  The parts std0 and std1 take about 2
# hours each, inv and mix 9 to 20 minutes, the others seconds.  Then run run_both_checks.sh.
cd "$(dirname "$0")"
mkdir -p logs
for fam in anti0 anti1 two0 two1 three0 three1 tto0 tto1 antistd0 antistd1 twostd0 twostd1 inv mix std0 std1; do
  bash guard.sh logs/both_$fam.log 36000 cert.py $fam 14+ 14+ save=logs/both_$fam.pkl
  echo "$fam: $(tail -n 1 logs/both_$fam.log) $(grep '^time' logs/both_$fam.log)"
done
