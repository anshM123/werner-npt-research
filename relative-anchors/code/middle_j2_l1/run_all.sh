#!/usr/bin/env bash
# Reruns every check and certificate of this folder, one Python process at a time, each under guard.sh (2.5 GB, time
# limit per step).  PASS or FAIL per step, ALL PASS at the end.  Steps marked "expect FAIL" are negative controls of
# certificates: they pass if the certificate prints FAIL (and the process was not stopped by the guard).
# Usage: bash run_all.sh > logs/run_all.out
cd "$(dirname "$0")"
mkdir -p logs
all=1
t_start=$(date +%s)
KMAX=${KMAX:-10}

step() {   # step NAME EXPECT LIMIT script args...
  local name="$1" expect="$2" lim="$3"; shift 3
  local log="logs/run_${name}.log"
  local t0=$(date +%s)
  bash guard.sh "$log" "$lim" "$@"
  local st=$?
  local last=$(tail -n 1 "$log" | tr -d '\r')
  local t1=$(date +%s)
  local res=FAIL
  if [ "$st" != 3 ]; then
    if [ "$expect" = PASS ] && [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; fi
    if [ "$expect" = FAIL ] && [ "$last" = FAIL ]; then res=PASS; fi
  fi
  [ "$res" = PASS ] || all=0
  echo "$res  $name  (python $*; expect $expect; $((t1 - t0)) s)"
}

echo "run_all.sh started $(date), KMAX = $KMAX"
# ---- certificates: every size pair (u0, u1), u0 >= u1, u0 + u1 = k, for 5 <= k <= KMAX (each saves what it certifies)
for k in $(seq 5 "$KMAX"); do
  for u1 in $(seq 0 $((k / 2))); do
    u0=$((k - u1))
    step "cert_k${k}_${u0}_${u1}" PASS 3600 cert_fixed.py "$u0" "$u1" "save=logs/cert_k${k}_${u0}_${u1}.pkl"
  done
done
# ---- coverage of all two-value weight vectors for 5 <= k <= KMAX
step coverage         PASS 60   coverage.py "$KMAX"
# ---- the builder, and the certified objects against j1lib.py (Lemma B) and ylib.py (definitions)
step builder          PASS 1800 check_builder.py
step ylib             PASS 1800 check_ylib3.py
step j1lib_k5_k6      PASS 1800 check_reduced.py 5,6 spectral=5,6
step j1lib_k7         PASS 1800 check_reduced.py 7 spectral=7
step j1lib_k8         PASS 3600 check_reduced.py 8
for k in $(seq 9 "$KMAX"); do
  step "j1lib_k${k}"  PASS 3600 check_reduced.py "$k"
done
# ---- negative controls: B -> (1 + 1/1000)B
for k in 5 6 7; do
  for u1 in $(seq 1 $((k / 2))); do
    u0=$((k - u1))
    step "neg_cert_k${k}_${u0}_${u1}" FAIL 1800 cert_fixed.py "$u0" "$u1" 1001/1000
  done
done
step neg_points       PASS 3600 neg_points.py
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
