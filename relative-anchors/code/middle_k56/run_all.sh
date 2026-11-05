#!/usr/bin/env bash
# Reruns every certificate and check of this folder (k = 5, then k = 6), one Python process at a time, each under guard.sh
# (2.5 GB, time limit per step).  PASS or FAIL per step, ALL PASS at the end.  Steps marked "expect FAIL" are negative
# controls of certificates: they pass if the certificate prints FAIL (and the process was not stopped by the guard).
# Usage: bash run_all.sh > logs/run_all.out
cd "$(dirname "$0")"
mkdir -p logs
all=1
t_start=$(date +%s)
FAMS="L3 T31 T32 T33 T42"

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

configs() {   # configs FAMILY [two]  -> "a,b" list (two: only configurations with two levels)
  python -c "
import k5lib as K
j, m = K.fam('$1')
print(' '.join(f'{a},{b}' for a, b in K.configs(5, m) if '$2' != 'two' or (a, b) not in ((0, 0), (m, 5 - m))))"
}

echo "run_all.sh started $(date)"
# ---- certificates: every family, every configuration (each saves the block and the minors it certifies)
for f in $FAMS; do
  for c in $(configs "$f"); do
    a=${c%,*}; b=${c#*,}
    step "cert_k5_${f}_${a}_${b}" PASS 1800 cert5.py "$f" "$a" "$b" "save=logs/cert_k5_${f}_${a}_${b}.pkl"
  done
done
# ---- coverage of every block and every two-value weight vector (exact identity at one ratio per pattern)
step coverage            PASS 3600 coverage.py L3,T31,T32,T33,T42
# ---- the certified blocks against j1lib.py (Lemma B) and ylib.py (definitions); minors; end to end
step builders_L3_T31     PASS 3600 check_builders.py L3,T31
step builders_T32_T33_T42 PASS 3600 check_builders.py T32,T33,T42
step minors              PASS 3600 check_minors.py L3,T31,T32,T33,T42
step end2end             PASS 3600 check_end2end.py L3,T31,T32,T33,T42
# ---- negative controls: B -> (1 + 1/1000)B
for f in $FAMS; do
  for c in $(configs "$f" two); do
    a=${c%,*}; b=${c#*,}
    step "neg_cert_k5_${f}_${a}_${b}" FAIL 1800 cert5.py "$f" "$a" "$b" 1001/1000
  done
done
step neg_points          PASS 3600 neg_points.py L3,T31,T32,T33,T42

# ======== k = 6: one certificate per configuration up to the level exchange (k5lib.rep_configs)
# direct certificates (cert5.py) for the smaller families, isotypic reduction (cert_sym.py) for the larger ones
D6="j4m4 j3m3 j5m3 j5m2 j4m3 j3m2 j4m2"
S6="j4m1 j4m0 j3m1 j3m0"
reps6() {   # reps6 FAMILY [two]
  python -c "
import k5lib as K
j, m = K.fam('$1')
print(' '.join(f'{a},{b}' for a, b in K.rep_configs(6, m) if '$2' != 'two' or (a, b) != (0, 0)))"
}
for f in $D6; do
  for c in $(reps6 "$f"); do
    a=${c%,*}; b=${c#*,}
    step "cert_k6_${f}_${a}_${b}" PASS 2400 cert5.py "$f" "$a" "$b" k=6 "save=logs/cert_k6_${f}_${a}_${b}.pkl"
  done
done
for f in $S6; do
  for c in $(reps6 "$f"); do
    a=${c%,*}; b=${c#*,}
    step "cert_k6_${f}_${a}_${b}" PASS 2400 cert_sym.py "$f" "$a" "$b" k=6 "save=logs/cert_k6_${f}_${a}_${b}.pkl"
  done
done
ALL6=$(echo $D6 $S6 | tr ' ' ',')
step coverage_k6         PASS 3600 coverage.py "$ALL6" k=6 sample=40
step builders_k6         PASS 3600 check_builders.py "$ALL6" k=6 reps=1
step minors_k6           PASS 3600 check_minors.py "$(echo $D6 | tr ' ' ',')" k=6 reps=1
step sym_k6              PASS 3600 check_sym.py "$(echo $S6 | tr ' ' ',')" k=6
step end2end_k6          PASS 3600 check_end2end.py "$ALL6" k=6 reps=1
for f in $D6; do
  for c in $(reps6 "$f" two); do
    a=${c%,*}; b=${c#*,}
    step "neg_cert_k6_${f}_${a}_${b}" FAIL 2400 cert5.py "$f" "$a" "$b" 1001/1000 k=6
  done
done
for f in $S6; do
  for c in $(reps6 "$f" two); do
    a=${c%,*}; b=${c#*,}
    step "neg_cert_k6_${f}_${a}_${b}" FAIL 2400 cert_sym.py "$f" "$a" "$b" 1001/1000 k=6
  done
done
step neg_points_k6       PASS 3600 neg_points.py "$ALL6" k=6 reps=1
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
