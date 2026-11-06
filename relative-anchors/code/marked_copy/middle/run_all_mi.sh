#!/usr/bin/env bash
# The new steps of NOTES_MI.md (Sections 2 to 6), one Python process at a time, each under guard.sh
# (2.5 GB, a time limit per step).  PASS or FAIL per step, ALL PASS at the end.  Steps marked "expect FAIL" are
# negative controls (B -> (1 + 1/1000) B): they pass if the certificate run ends with FAIL (and the process was not
# stopped by the guard).  The older results used by the coverage (Theorems M, B, F, I, N1 to N4 and the fixed-k
# certificates for 33 <= k <= 40) are read from logs/*.pkl and logs/frozen/, written by run_all.sh, run_extra.sh,
# run_extra3.sh and run_extra4.sh of this folder (run them first).
# Usage: bash run_all_mi.sh > logs/run_all_mi.out
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

echo "run_all_mi.sh started $(date)"
export GAUGE=R
OLD="logs/mid_i1.pkl logs/mid_i2.pkl logs/mid_i3.pkl logs/mid_i4.pkl logs/fam_1_1.pkl logs/fam_lowirr_0.pkl logs/frozen/band_bottom_8.pkl logs/frozen/band_top_8.pkl"
EDGEQ="-1,n,n 0,n,n 1,n,n 2,n,n 3,n,n"
EDGET="n,-2,n n,-1,n n,0,n n,1,n n,2,n n,3,n"
# ---- the generic family [q0 >= 4, t >= 4, i >= 5]: the (u - i - 1, i, 1) part (four variables x, y, z, sigma)
step gen3             PASS 1800 prove_gen.py logs/gen3.pkl n,n,n --imin=5 --parts=3
# ---- the edge families with i symbolic, whole form: [q0 = a, t >= 4, i >= 6], [q0 >= 4, t = c, i >= 5]
step edges_q          PASS 5400 prove_gen.py logs/edges_q.pkl $EDGEQ --imin=6
step edges_t          PASS 5400 prove_gen.py logs/edges_t.pkl $EDGET --imin=5
# ---- the family [q0 >= 4, t >= 4, i = 5], whole form (as Theorems N1 to N4)
step mid_i5           PASS 3000 prove_gen.py logs/mid_i5.pkl n,n,5
# ---- fixed-k certificates for the forms with 41 <= k <= 50 outside the proved families
NEW="logs/gen3.pkl logs/edges_q.pkl logs/edges_t.pkl logs/mid_i5.pkl"
step cert_rest_mi     PASS 3600 cert_rest_mi.py 41 50 logs/members_41_50.pkl $NEW $OLD
# ---- coverage of every (k, j, i) with k <= 120
step coverage_mi      PASS 3600 coverage_mi.py 120 32 6 51 logs/members_33_40.pkl logs/members_41_50.pkl $NEW $OLD
# ---- the new family forms against mlib.py at far members (k up to about 80)
step famcheck_gen     PASS 1800 famcheck_mi.py n,n,n n,n,5 $EDGET --imin=5
step famcheck_edgeq   PASS 1800 famcheck_mi.py $EDGEQ --imin=6
# ---- exact eliminations at members with 41 <= k <= 60 (and negative controls there)
step memcheck         PASS 3600 memcheck_mi.py 41 60 7 n,n,n $EDGEQ $EDGET --imin=5
# ---- negative controls of the certificates, B -> (1 + 1/1000) B
step neg_gen3         FAIL 1800 prove_gen.py logs/unused.pkl n,n,n --imin=5 --parts=3 --bscale=1001/1000
step neg_edges_q      FAIL 3600 prove_gen.py logs/unused.pkl $EDGEQ --imin=6 --bscale=1001/1000 --tlimit=300
step neg_edges_t      FAIL 3600 prove_gen.py logs/unused.pkl $EDGET --imin=5 --bscale=1001/1000 --tlimit=300
step neg_i5           FAIL 1800 prove_gen.py logs/unused.pkl n,n,5 --bscale=1001/1000 --tlimit=300
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
