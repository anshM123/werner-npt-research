#!/usr/bin/env bash
# Extension (run it after run_all_mi.sh, which produces the pickles used here): the family [q0 >= 4, t >= 4, i = 6],
# fixed-k certificates for 51 <= k <= 55, the coverage with them, a far-member check and the negative control.  One
# Python process at a time, each under guard.sh; PASS or FAIL per step, ALL PASS at the end.
# Usage: bash run_extra_mi.sh > logs/run_extra_mi.out
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

echo "run_extra_mi.sh started $(date)"
export GAUGE=R
OLD="logs/mid_i1.pkl logs/mid_i2.pkl logs/mid_i3.pkl logs/mid_i4.pkl logs/fam_1_1.pkl logs/fam_lowirr_0.pkl logs/frozen/band_bottom_8.pkl logs/frozen/band_top_8.pkl"
NEW="logs/gen3.pkl logs/edges_q.pkl logs/edges_t.pkl logs/mid_i5.pkl"
# ---- the family [q0 >= 4, t >= 4, i = 6], whole form
step mid_i6           PASS 3000 prove_gen.py logs/mid_i6.pkl n,n,6
# ---- fixed-k certificates for the forms with 51 <= k <= 55 outside the proved families
step cert_rest_mi2    PASS 2400 cert_rest_mi.py 51 55 logs/members_51_55.pkl logs/mid_i6.pkl $NEW $OLD
# ---- coverage of every (k, j, i) with k <= 120: what is not covered has i >= 7 and k >= 56
step coverage_mi2     PASS 3600 coverage_mi.py 120 32 7 56 logs/members_33_40.pkl logs/members_41_50.pkl logs/members_51_55.pkl logs/mid_i6.pkl $NEW $OLD
# ---- the family form against mlib.py at far members
step famcheck_i6      PASS 1800 famcheck_mi.py n,n,6
# ---- negative control, B -> (1 + 1/1000) B
step neg_i6           FAIL 1800 prove_gen.py logs/unused.pkl n,n,6 --bscale=1001/1000 --tlimit=300
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
