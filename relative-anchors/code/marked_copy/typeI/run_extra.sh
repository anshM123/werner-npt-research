#!/usr/bin/env bash
# The extension (Section 8 of NOTES.md), run after run_all.sh: the fixed-k certificates for 21 <= k <= 28, the
# coverage for 8 <= k <= 28 (with the files of run_all.sh) and the bookkeeping of every block type, the certificates for
# k = 4..7 (also covered by Theorems 5_2, 6_2, 7_2), symbolic against numeric forms at every 20th form for
# 18 <= k <= 24, and whole blocks end to end at k = 10 and 12 for some j.  One Python process at a time under
# guard.sh; PASS or FAIL per step, ALL PASS at the end.
# Usage: bash run_extra.sh > logs/run_extra.out
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
  local last=$(tail -n 1 "$log" | tr -d '\r' | awk '{print $1}')
  local t1=$(date +%s)
  local res=FAIL
  if [ "$st" != 3 ]; then
    if [ "$expect" = PASS ] && [ "$st" = 0 ] && [ "$last" = PASS ]; then res=PASS; fi
    if [ "$expect" = FAIL ] && [ "$last" = FAIL ]; then res=PASS; fi
  fi
  [ "$res" = PASS ] || all=0
  echo "$res  $name  (python $*; expect $expect; $((t1 - t0)) s)"
}

echo "run_extra.sh started $(date)"
step cert_04_07       PASS 1800 cert_t.py 4 7 1 logs/cert_04_07.pkl
step coverage_4_7     PASS 1800 coverage_t.py 4 7 logs/cert_04_07.pkl
step symnum_far       PASS 1800 check_symnum.py 18 24 20
step e2e_10           PASS 2400 e2e_t.py 10 1/2 3,8,9
step e2e_12           PASS 3000 e2e_t.py 12 2 3,10,11
for k in 21 22 23 24 25 26 27 28; do
  step cert_$k        PASS 3600 cert_t.py $k $k 1 logs/cert_$k.pkl
done
PK="logs/cert_08_12.pkl logs/cert_13_15.pkl logs/cert_16_17.pkl logs/cert_18.pkl logs/cert_19.pkl logs/cert_20.pkl"
for k in 21 22 23 24 25 26 27 28; do PK="$PK logs/cert_$k.pkl"; done
step coverage_8_28    PASS 1800 coverage_t.py 8 28 $PK
step blocks_8_28      PASS 1800 blocks_t.py 8 8 28 $PK
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
