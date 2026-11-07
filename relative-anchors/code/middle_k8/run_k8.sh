#!/usr/bin/env bash
# k = 8 at two weight values: the fast Lemma B builder against j1lib, the certificates for the 26 open families
# (reduced method), then the checks and the second check.  One Python process at a time under guard.sh.
# Usage: bash run_k8.sh > logs/run_k8.out
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
FAMS8="j3m0 j3m1 j3m2 j3m3 j4m0 j4m1 j4m2 j4m3 j4m4 j5m0 j5m1 j5m2 j5m3 j5m4 j5m5 j6m0 j6m1 j6m2 j6m3 j6m4 j6m5 j6m6 j7m2 j7m3 j7m4 j7m5"
reps8() {   # reps8 FAMILY [two]
  python -c "
import k5lib as K
j, m = K.fam('$1')
print(' '.join(f'{a},{b}' for a, b in K.rep_configs(8, m) if '$2' != 'two' or (a, b) != (0, 0)))"
}
echo "run_k8.sh started $(date)"
# the fast Lemma B builder (used by every check below through k5lib.build_j1) against j1lib.J1, every torus group
step fastpd    PASS 1800 validate_pd.py 30 11
step fastj1_k4 PASS 1800 validate.py 4 2 2
step fastj1_k5 PASS 3600 validate.py 5 2 3
step fastj1_k6 PASS 7200 validate.py 6 1 4
for f in $FAMS8; do
  for c in $(reps8 "$f"); do
    a=${c%,*}; b=${c#*,}
    step "cert_k8_${f}_${a}_${b}" PASS 7200 cert_sym.py "$f" "$a" "$b" k=8 "save=logs/cert_k8_${f}_${a}_${b}.pkl"
  done
done
ALL8=$(echo $FAMS8 | tr ' ' ',')
step coverage_k8   PASS 21600 coverage.py "$ALL8" k=8 sample=200
step sym_k8        PASS 21600 check_sym.py "$ALL8" k=8
step builders_k8   PASS 21600 check_builders.py "$ALL8" k=8 reps=1
for f in $FAMS8; do
  for c in $(reps8 "$f" two); do
    a=${c%,*}; b=${c#*,}
    step "neg_cert_k8_${f}_${a}_${b}" FAIL 7200 cert_sym.py "$f" "$a" "$b" 1001/1000 k=8
  done
done
step neg_points_k8 PASS 21600 neg_points.py "$ALL8" k=8 reps=1
step end2end_k8    PASS 21600 check_end2end.py "$ALL8" k=8 reps=1
# the second check: every torus block of every (R_inf)_j at random two-value points; the saved blocks at 16 values of q
step crosscheck_neg      PASS 900   crosscheck.py N
step crosscheck_blocks   PASS 21600 crosscheck.py A 8 12 8208
step crosscheck_saved    PASS 21600 crosscheck.py C logs 31
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
