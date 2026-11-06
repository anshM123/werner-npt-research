#!/usr/bin/env bash
# Stage 1 extended: the regimes (b+, b) for 31 <= b <= 40 (the larger level symbolic, the smaller fixed at b),
# all 16 families each (batch.py), one Python process at a time under guard.sh.
cd "$(dirname "$0")"
mkdir -p logs
for b in $(seq 31 40); do
  bash guard.sh logs/ext_cert_${b}.log 3600 batch.py "${b}+,${b}"
  echo "cert_${b}+_${b}: $(tail -n 1 logs/ext_cert_${b}.log) $(grep -c 'PASS' logs/ext_cert_${b}.log) lines PASS"
done
