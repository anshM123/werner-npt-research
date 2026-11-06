#!/usr/bin/env bash
# Reruns every check and certificate of NOTES.md, one Python process at a time under guard.sh (2.5 GB limit).
# Usage: bash run_all.sh > logs/run_all.out     (M = 13: the stage-1 regimes (o(m)+, m), 2 <= m <= 13)
cd "$(dirname "$0")"
mkdir -p logs
M=${M:-13}
allok=1
t00=$(date +%s)
step() {   # step NAME EXPECT TIMEOUT args...
  name="$1"; shift; expect="$1"; shift; lim="$1"; shift
  t0=$(date +%s)
  bash guard.sh "logs/run_${name}.log" "$lim" "$@"
  st=$?
  last=$(tail -n 1 "logs/run_${name}.log" | tr -d '\r')
  got=FAIL
  if [ $st -eq 0 ] && [ "$last" = "PASS" ]; then got=PASS; fi
  res=PASS
  if [ "$got" != "$expect" ]; then res=FAIL; allok=0; fi
  echo "$res  $name  (python $*; expect $expect; $(( $(date +%s) - t0 )) s)"
}
# 1. the reduction against the explicit block of the Lemma B builder (k = 5, ..., 8, and k = 11 inside the regimes)
step points_k5 PASS 600 check_points.py 3 2 1/9 4
step points_k5b PASS 600 check_points.py 2 3 49/25
step points_k6 PASS 600 check_points.py 4 2 1/4 9/4
step points_k6b PASS 600 check_points.py 3 3 1/9
step points_k6c PASS 600 check_points.py 2 4 4
step points_k7 PASS 900 check_points.py 5 2 4
step points_k7b PASS 900 check_points.py 4 3 1/9
step points_k7c PASS 900 check_points.py 3 4 9/4
step points_k8 PASS 900 check_points.py 6 2 1/4
step points_k8b PASS 900 check_points.py 5 3 49/25
step points_k8c PASS 900 check_points.py 4 4 1/9
step points_k11 PASS 1800 check_points.py 9 2 1/9 sample=12,2
step points_k11b PASS 1800 check_points.py 7 4 49/25 sample=10,4
step points_k11c PASS 1800 check_points.py 6 5 4 sample=9,5
step points_k12 PASS 1800 check_points.py 8 4 1/4 sample=10,4
# 2. the certificates (stage 1: one level of fixed size m, the other symbolic)
regs=""
for m in $(seq 2 "$M"); do
  o=$(( m > 11 - m ? m : 11 - m ))
  regs="$regs ${o}+,${m}"
done
for r in $regs; do
  step "cert_${r/,/_}" PASS 1800 batch.py "$r"
done
# 3. the certified forms against the rational-arithmetic reduction at points
step forms PASS 1800 check_forms.py 9+,2 8+,3 7+,4 6+,5 "${M}+,${M}"
# 4. negative controls
step neg PASS 900 neg.py
# 5. coverage of all two-value configurations with 11 <= k <= 60
step coverage PASS 900 coverage.py "$M" none 60
step coverage_k80 PASS 900 coverage.py "$M" none 80
echo "total $(( $(date +%s) - t00 )) s"
if [ $allok = 1 ]; then echo "ALL PASS"; else echo "SOME FAIL"; fi
