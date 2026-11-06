#!/usr/bin/env bash
# k = 7 at two weight values: certificates for the 18 open families (reduced method), then the checks.
# One Python process at a time under guard.sh.  Usage: bash run_k7.sh > logs/run_k7.out
cd "$(dirname "$0")"
mkdir -p logs
all=1
t_start=$(date +%s)
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
FAMS7="j3m0 j3m1 j3m2 j3m3 j4m0 j4m1 j4m2 j4m3 j4m4 j5m0 j5m1 j5m2 j5m3 j5m4 j5m5 j6m2 j6m3 j6m4"
reps7() {   # reps7 FAMILY [two]
  python -c "
import k5lib as K
j, m = K.fam('$1')
print(' '.join(f'{a},{b}' for a, b in K.rep_configs(7, m) if '$2' != 'two' or (a, b) != (0, 0)))"
}
echo "run_k7.sh started $(date)"
for f in $FAMS7; do
  for c in $(reps7 "$f"); do
    a=${c%,*}; b=${c#*,}
    step "cert_k7_${f}_${a}_${b}" PASS 3600 cert_sym.py "$f" "$a" "$b" k=7 "save=logs/cert_k7_${f}_${a}_${b}.pkl"
  done
done
ALL7=$(echo $FAMS7 | tr ' ' ',')
step coverage_k7   PASS 14400 coverage.py "$ALL7" k=7 sample=800
step sym_k7        PASS 14400 check_sym.py "$ALL7" k=7
step builders_k7   PASS 14400 check_builders.py "$ALL7" k=7 reps=1
for f in $FAMS7; do
  for c in $(reps7 "$f" two); do
    a=${c%,*}; b=${c#*,}
    step "neg_cert_k7_${f}_${a}_${b}" FAIL 3600 cert_sym.py "$f" "$a" "$b" 1001/1000 k=7
  done
done
step neg_points_k7 PASS 14400 neg_points.py "$ALL7" k=7 reps=1
step end2end_k7    PASS 21600 check_end2end.py "$ALL7" k=7 reps=1
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
