#!/usr/bin/env bash
# Second extension (run it after run_all.sh and run_extra.sh, which produce the pickles used here): the family
# [q0 >= 4, t >= 4, i = 3], the coverage with it, and its negative controls.  One Python process at a time, each under
# guard.sh; PASS or FAIL per step, ALL PASS at the end.
# Usage: bash run_extra3.sh > logs/run_extra3.out
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

echo "run_extra3.sh started $(date)"
export GAUGE=R
FAMS="logs/mid_i1.pkl logs/fam_1_1.pkl logs/fam_lowirr_0.pkl logs/frozen/band_bottom_8.pkl logs/frozen/band_top_8.pkl"
# ---- the family [q0 >= 4, t >= 4, i = 3]
step mid_i3           PASS 2400 prove_mid.py logs/mid_i3.pkl n,n,3
# ---- coverage with it: what is not covered lies in 9 <= j <= k - 9, i >= 4, k >= 41
step coverage_i3      PASS 1800 coverage_mid.py 80 32 band:4 logs/members_33_40.pkl logs/mid_i3.pkl logs/mid_i2.pkl $FAMS
# ---- negative controls, B -> (1 + 1/1000) B
step neg_i3           FAIL 1800 prove_mid.py logs/unused.pkl n,n,3 --bscale=1001/1000
step neg_points_i3    PASS 1800 neg_mid.py 40 n,n,3
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
