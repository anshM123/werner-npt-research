#!/usr/bin/env bash
# development test of the k = 6 check steps on j4m4 (direct) and j4m1 (reduced), from the development certificates
cd "$(dirname "$0")/../.."
for f in j4m4 j4m1; do
  for p in logs/dev/cert_k6_${f}_*.pkl; do cp "$p" logs/; done
  for l in logs/dev/cert_k6_${f}_*.log; do b=$(basename "$l"); cp "$l" "logs/run_${b}"; done
done
for cmd in "coverage.py j4m4,j4m1 k=6 sample=40" "check_builders.py j4m4,j4m1 k=6 reps=1" "check_minors.py j4m4 k=6 reps=1" \
           "check_sym.py j4m1 k=6" "check_end2end.py j4m4,j4m1 k=6 reps=1" "neg_points.py j4m4,j4m1 k=6 reps=1" \
           "cert_sym.py j4m1 0 2 1001/1000 k=6"; do
  name=$(echo "$cmd" | cut -d' ' -f1 | sed 's/.py//')
  bash guard.sh "logs/dev/t6_${name}.log" 1800 $cmd
  echo "== $cmd -> exit $?: $(tail -n1 logs/dev/t6_${name}.log | tr -d '\r')"
  grep -E "^[a-zA-Z0-9]+:|MISMATCH|FAIL|NOT|Error|Traceback" "logs/dev/t6_${name}.log" | head -8
done
echo TESTS DONE
