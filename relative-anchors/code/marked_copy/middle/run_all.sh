#!/usr/bin/env bash
# Reruns every check and certificate of this folder, one Python process at a time, each under guard.sh (2.5 GB, time
# limit per step).  PASS or FAIL per step, ALL PASS at the end.  Steps marked "expect FAIL" are negative controls of
# certificates: they pass if the certificate run ends with FAIL (and the process was not stopped by the guard).
# Theorem M (fixed k <= 32) and the band families of Theorem B are results of marked_copy; their saved outcomes are read
# from logs/frozen/ (copied unchanged from ../logs, whose run_all ended ALL PASS).  Theorem F (one symbolic
# size) and the family of Theorem I are rerun here (they are cheap) and also rerun with the split certificate.
# Usage: bash run_all.sh > logs/run_all.out
cd "$(dirname "$0")"
mkdir -p logs
all=1
t_start=$(date +%s)
ONESIZE="-1,-2,n -1,-1,n -1,0,n -1,1,n -1,2,n -1,3,n -1,n,0 -1,n,1 0,-2,n 0,-1,n 0,0,n 0,1,n 0,2,n 0,3,n 0,n,0 0,n,1 1,-2,n 1,-1,n 1,0,n 1,1,n 1,2,n 1,3,n 1,n,0 1,n,1 2,-2,n 2,-1,n 2,0,n 2,1,n 2,2,n 2,3,n 2,n,0 2,n,1 3,-2,n 3,-1,n 3,0,n 3,1,n 3,2,n 3,3,n 3,n,0 3,n,1 n,-2,0 n,-2,1 n,-1,0 n,-1,1 n,0,0 n,0,1 n,1,0 n,1,1 n,2,0 n,2,1 n,3,0 n,3,1"

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

echo "run_all.sh started $(date)"
export GAUGE=R
# ---- the orbit calculus against the Lemma B builder (check of marked_copy): k <= 6 exact; symbolic against numeric
step builder          PASS 1800 check_builder.py 6 6 16
# ---- the irrep split against the Lemma B builder, class sums, kernel vectors (k <= 7)
step split            PASS 1800 check_split.py 7
# ---- the families of marked_copy used by the coverage: one symbolic size (Theorem F), [q0>=4, t>=4, i=0] (Theorem I)
step fam1             PASS 3600 prove_fam.py 1 1
step fam_invariant    PASS 1800 prove_lowirr.py 0
# ---- the split certificate on the same 52 families (a second, independent elimination)
step fam1_split       PASS 3600 prove_mid.py logs/fam1_split.pkl $ONESIZE
# ---- new: the family [q0 >= 4, t >= 4, i = 1]
step mid_i1           PASS 5400 prove_mid.py logs/mid_i1.pkl n,n,1
# ---- fixed-k certificates (cert_fixed.certify of marked_copy) for every form with 33 <= k <= 40 outside the families
FAMS="logs/mid_i1.pkl logs/fam_1_1.pkl logs/fam_lowirr_0.pkl logs/frozen/band_bottom_8.pkl logs/frozen/band_top_8.pkl"
step cert_rest        PASS 3600 cert_rest.py 33 40 logs/members_33_40.pkl $FAMS
# ---- coverage of every (k, j, i) with k <= 80: what is not covered lies in 9 <= j <= k - 9, i >= 2, k >= 41
step coverage         PASS 1800 coverage_mid.py 80 32 band:2 logs/members_33_40.pkl $FAMS
# ---- negative controls: B -> (1 + 1/1000) B; the split certificate must fail on the one-size families and on the
#      new family; explicit members with 33 <= k <= 40 of the new family where the scaled form is not positive definite
step neg_split        FAIL 3600 prove_mid.py logs/unused.pkl $ONESIZE --bscale=1001/1000
step neg_i1           FAIL 1800 prove_mid.py logs/unused.pkl n,n,1 --bscale=1001/1000
step neg_points       PASS 1800 neg_mid.py 40 n,n,1
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
