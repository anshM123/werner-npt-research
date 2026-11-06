#!/usr/bin/env bash
# The second check of the marked-copy theorems (markcheck.py, projcheck.py, famcheck.py), one Python process at a
# time under guard.sh.  Usage: bash run_check.sh > logs/run_check.out
cd "$(dirname "$0")"
mkdir -p logs
all=1
t_start=$(date +%s)
step() {   # step NAME LIMIT script args...
  local name="$1" lim="$2"; shift 2
  local log="logs/run_check_${name}.log"
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
echo "run_check.sh started $(date)"
# the fast builder against Lemma B (J1.matrix), entry by entry
step builder_k4      600  markcheck.py V 4 1/3,2
step builder_k5      600  markcheck.py V 5 2/5,3
step builder_k6      900  markcheck.py V 6 1/3,5/2
# negative controls: B scaled by 1001/1000 must be refused somewhere (and, at k = 5, equals J1 with bscale)
step negative_k5     600  markcheck.py N 5 1/1000,1/3
step negative_k7     600  markcheck.py N 7 1/1000
# end to end, exact: every j, several sigma
step blocks_k7       1800 markcheck.py X 7 1/10,1/3,7/5,2,20
step blocks_k8       3000 markcheck.py X 8 1/10,1/3,7/5,2,20
step blocks_k9       3000 markcheck.py X 9 1/10,2,20
# the reduced forms of mlib.py against exact projections of the Lemma B block
step reduced_k7      1800 projcheck.py 7 2,1/3,1/10
step reduced_k8      2400 projcheck.py 8 2,1/3,1/10
step reduced_k9      3000 projcheck.py 9 2,1/10
# the family forms of famlib.py against mlib.py at members far from the sample points
step families_inv    1200 famcheck.py inv
step families_one    3600 famcheck.py fam1
step band_bottom     3600 famcheck.py bottom 8
step band_top        3600 famcheck.py top 8
# the two-size families of the middle band (marked_copy/middle) against mlib.py at far members, and their
# conclusion at members with 41 <= k <= 60 (exact rank and positive definiteness of a complement)
step families_mid    3600 famcheck.py mid
step members_41_60   1800 mbcheck.py 40 4141
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
