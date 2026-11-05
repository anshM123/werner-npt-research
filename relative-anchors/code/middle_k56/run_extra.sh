#!/usr/bin/env bash
# A second, independent certificate (not needed for any proof): the isotypic reduction of cert_sym.py, with the
# complement W in place of the deletion gauge, for every configuration that run_all.sh certifies directly
# (k = 5: all 52; k = 6: the 56 representative configurations of the direct families); and denser samples of the
# exact identity of coverage.py at k = 6.  One Python process at a time under guard.sh; PASS or FAIL per step, ALL PASS
# at the end.
# Usage: bash run_extra.sh > logs/run_extra.out
cd "$(dirname "$0")"
mkdir -p logs/extra
all=1
t_start=$(date +%s)

step() {   # step NAME LIMIT script args...
  local name="$1" lim="$2"; shift 2
  local log="logs/extra/run_${name}.log"
  local t0=$(date +%s)
  bash guard.sh "$log" "$lim" "$@"
  local st=$?
  local last=$(tail -n 1 "$log" | tr -d '\r')
  local t1=$(date +%s)
  local res=FAIL
  if [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; fi
  [ "$res" = PASS ] || all=0
  echo "$res  $name  (python $*; $((t1 - t0)) s)"
}

echo "run_extra.sh started $(date)"
for f in L3 T31 T32 T33 T42; do
  for c in $(python -c "
import k5lib as K
j, m = K.fam('$f')
print(' '.join(f'{a},{b}' for a, b in K.configs(5, m)))"); do
    a=${c%,*}; b=${c#*,}
    step "sym_k5_${f}_${a}_${b}" 1800 cert_sym.py "$f" "$a" "$b" k=5
  done
done
for f in j4m4 j3m3 j5m3 j5m2 j4m3 j3m2 j4m2; do
  for c in $(python -c "
import k5lib as K
j, m = K.fam('$f')
print(' '.join(f'{a},{b}' for a, b in K.rep_configs(6, m)))"); do
    a=${c%,*}; b=${c#*,}
    step "sym_k6_${f}_${a}_${b}" 1800 cert_sym.py "$f" "$a" "$b" k=6
  done
done
# denser samples of the exact coverage identity at k = 6 (run_all.sh checks every 40th pattern)
step coverage_k6_sample8 3600 coverage.py j4m4,j3m3,j5m3,j5m2,j4m3,j3m2,j4m2,j4m1,j3m1 k=6 sample=8
step coverage_k6_letter1 3600 coverage.py j4m0,j3m0 k=6 sample=4
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
