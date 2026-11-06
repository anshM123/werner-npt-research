#!/usr/bin/env bash
# The letter-1 block of (R_inf)_2 at two weight values with BOTH level sizes symbolic: u0 = 14 + x0, u1 = 14 + x1
# (x0, x1 >= 0 real), certificates in (q, x0, x1) for every family.  Two lanes, one Python process per lane.
# Usage: bash run_both.sh A|B
cd "$(dirname "$0")"
mkdir -p logs
if [ "$1" = A ]; then fams="anti0 anti1 two0 two1 three0 three1 tto0 tto1 antistd0 antistd1 twostd0 twostd1 std0 inv"; else fams="mix std1"; fi
for fam in $fams; do
  bash guard.sh logs/both_$fam.log 21600 cert.py $fam 14+ 14+ save=logs/both_$fam.pkl
  echo "$fam: $(tail -n 1 logs/both_$fam.log) $(grep '^time' logs/both_$fam.log)"
done
