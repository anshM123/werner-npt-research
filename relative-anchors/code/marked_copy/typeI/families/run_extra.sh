#!/usr/bin/env bash
# Run after run_all.sh (NOTES.md, Sections 1 and 4 to 7): the check behind the anti reduction to i = 2, the
# determinant route for the sym families of case MT with two symbolic sizes (Proposition D for i and tau'
# symbolic, Proposition D1 with slices in one size for the others), the measurements for the open families, and
# the coverage with every result.
# One Python process at a time, each under guard.sh.  Usage: bash run_extra.sh > logs/run_extra.out
cd "$(dirname "$0")"
mkdir -p logs res
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
# the two-symbolic MT sym families other than (i, tau'): with i = 0 or tau' = 0 (small) and the others (large)
tr -d '\r' < lists/MT_sym_2e.txt | awk -F'[/,]' '($5 == "0" || $6 == "0")' > lists/MT_sym_2e_small.txt
tr -d '\r' < lists/MT_sym_2e.txt | awk -F'[/,]' '!($5 == "0" || $6 == "0")' > lists/MT_sym_2e_large.txt
# ---- the anti reduction to i = 2 (Section 1): for every anti family with i and tau' symbolic, its form and every
#      intermediate denominator of its calculus depend on (i, tau') only through nu
step structure_anti PASS 1800 structure.py --dens @lists/MT_anti_all.txt @lists/MU_anti_all.txt
# ---- the determinant route: MT sym [q0 = a, t = c, i >= 2, tau' >= 2], every (a, c) of the list
step det_MT         PASS 5400 tdet.py prove res/det_MT.pkl MT $(cat lists/MT_sym_2d.txt | tr -d '\r')
# ---- measurements of the determinant route for families with three and four symbolic sizes
step measure_MT     PASS 1800 tdet.py measure MT n n
step measure_MT_q   PASS 1800 tdet.py measure MT 1 n
step measure_MT_t   PASS 1800 tdet.py measure MT n 1
step measure_MU     PASS 3600 tdet.py measure MU n n
step measure_p1     PASS 1800 tdet.py measure_plain MT/sym/n,n,1,0
step measure_p2     PASS 1800 tdet.py measure_plain MT/sym/1,n,n,0
step measure_p3     PASS 1800 tdet.py measure_plain MU/sym/n,n,0,1
# ---- the two-symbolic MT sym families with i = 0 or tau' = 0 (at most 10 vectors): the determinant route with
#      slices in one size (Proposition D1); the other 19 (13 vectors) are left open (NOTES.md, Section 6)
step det1_MT        PASS 3600 tdet1.py res/det1_MT.pkl $(cat lists/MT_sym_2e_small.txt | tr -d '\r')
# ---- coverage with every result
step coverage_all   PASS 3600 coverage_st3.py 120 29 30 res/anti_MT.pkl res/anti_MU.pkl res/sym1_MT.pkl res/sym1_MU.pkl res/det_MT.pkl res/det1_MT.pkl
# ---- summary of the certified families
step summary        PASS 600  summary.py res/anti_MT.pkl res/anti_MU.pkl res/sym1_MT.pkl res/sym1_MU.pkl res/det_MT.pkl res/det1_MT.pkl
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
