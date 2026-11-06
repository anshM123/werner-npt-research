#!/usr/bin/env bash
# The steps of NOTES.md, one Python process at a time, each under guard.sh (2.0 GB, a time limit per step).
# PASS or FAIL per step, ALL PASS at the end.  Steps marked "expect FAIL" are negative controls (B -> (1 + 1/1000) B):
# they pass if the run ends with FAIL (and the process was not stopped by the guard).
# Usage: bash run_all.sh > logs/run_all.out
cd "$(dirname "$0")"
mkdir -p logs res lists
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

echo "run_all.sh started $(date)"
# ---- the family lists
step lists          PASS 900  mklists.py
# ---- the structure lemma (dependence on (i, tau') through nu = tau' + 2i) on generic and edge forms
step structure      PASS 3600 structure.py MT/anti/n,n,n,n MU/anti/n,n,n,n MT/anti/0,1,n,n MU/anti/1,-1,n,n MT/sym/0,1,n,n MU/sym/1,-1,n,n MT/sym/2,n,n,n MU/sym/n,2,n,n MT/sym/n,n,n,n MU/sym/n,n,n,n
# ---- the family forms against tlib.py at members far from the sanity points (k <= 70)
step famcheck       PASS 3600 tfamcheck.py 70 MT/anti/n,n,n,n MU/anti/n,n,n,n MT/anti/n,n,1,n MT/anti/1,n,0,n MU/anti/n,-2,n,n MU/anti/2,n,1,n MT/sym/n,1,1,0 MT/sym/2,n,0,1 MT/sym/1,2,n,0 MT/sym/2,3,1,n MT/sym/0,1,n,n MT/sym/n,n,1,0 MT/sym/2,n,n,1 MT/sym/n,2,1,n MU/sym/n,1,0,1 MU/sym/2,n,1,1 MU/sym/-1,2,n,1 MU/sym/1,-1,1,n MU/sym/1,0,n,n MU/sym/n,n,1,1 MU/sym/0,n,n,1 MU/sym/n,-2,n,1 MT/sym/n,n,n,n MU/sym/n,n,n,n
# ---- certificates: every anti family (all symbolic-size patterns)
step anti_MT        PASS 1800 tprove.py res/anti_MT.pkl @lists/MT_anti_all.txt
step anti_MU        PASS 3600 tprove.py res/anti_MU.pkl @lists/MU_anti_all.txt
# ---- certificates: the sym families with one symbolic size
step sym1_MT        PASS 1800 tprove.py res/sym1_MT.pkl @lists/MT_sym_1.txt
step sym1_MU        PASS 5400 tprove.py res/sym1_MU.pkl @lists/MU_sym_1.txt
# ---- coverage for k <= 120 (Theorem I1 for k <= 29, the families above)
step coverage       PASS 3600 coverage_st3.py 120 29 30 res/anti_MT.pkl res/anti_MU.pkl res/sym1_MT.pkl res/sym1_MU.pkl
# ---- negative controls: B -> (1 + 1/1000) B (expect FAIL: some certificate fails)
step neg_sym1_MT    FAIL 1800 tprove.py - @lists/MT_sym_1.txt --bscale=1001/1000
step neg_anti       FAIL 1800 tprove.py - @lists/MT_anti_4.txt @lists/MU_anti_4.txt --bscale=1001/1000
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
