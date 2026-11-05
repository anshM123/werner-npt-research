#!/usr/bin/env bash
# Reruns every check and certificate of this folder, one Python process at a time, each under guard.sh (2.5 GB, time
# limit per step).  PASS or FAIL per step, ALL PASS at the end.  Steps marked "expect FAIL" are negative controls of
# certificates: they pass if the certificate prints FAIL (and the process was not stopped by the guard).
# The two regimes whose certificates did not finish (NOTES.md, rows 5b and 10b) are not run here; see attempts.sh.
# Usage: bash run_all.sh > logs/run_all.out
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

M2=12   # family 2: the inv part is certified at the fixed-size regimes (m, 3+) and (3+, m), 3 <= m <= M2
M1=12   # family 1: the regime tw = (0, 1) is certified at the fixed-size regimes (m, 2+), 2 <= m <= M1

echo "run_all.sh started $(date)"
# ---- block structure and builders
step kernels         PASS 900  explore_kernel.py 6
step kernels_k7      PASS 1200 explore_kernel.py 7
step ylib            PASS 1800 check_ylib_r2.py
step ylib_reduced    PASS 900  check_ylib_reduced.py
# ---- family 1: |T| = 2
step f1_check        PASS 600  check_t2.py
step f1_cert_one     PASS 600  cert_t2.py 2+
step f1_cert_A_1_2   PASS 1800 cert_t2.py 1 2+ tw=0,0
step f1_cert_A_2_1   PASS 1800 cert_t2.py 2+ 1 tw=0,0
step f1_cert_A_2_2   PASS 3600 cert_t2.py 2+ 2+ tw=0,0
step f1_cert_B_1_2   PASS 1800 cert_t2.py 1 2+ tw=0,1
for m in $(seq 2 "$M1"); do
  step f1_cert_B_${m}_2 PASS 900 cert_t2.py "$m" 2+ tw=0,1
done
step f1_neg_points   PASS 1800 neg_t2.py
step f1_neg_cert_one FAIL 600  cert_t2.py 2+ form=G bscale=1001/1000
step f1_gauge_cert_one PASS 600 cert_t2.py 2+ form=G
# ---- family 2: |T| = 1
step f2_check        PASS 900  check_f2.py
step f2_dims         PASS 600  dims.py
step f2_one_inv      PASS 600  cert_f2.py inv 3+ 0
step f2_one_std      PASS 600  cert_f2.py std 3+ 0
step f2_one_two      PASS 600  cert_f2.py two 4+ 0
step f2_one_anti     PASS 600  cert_f2.py anti 3+ 0
step f2_neg_cert_one FAIL 600  cert_f2.py inv 3+ 0 bscale=1001/1000
for tw in 0 1; do
  step f2_std_2_2_tw$tw    PASS 1800 cert_f2.py std 2 2 tw=$tw
  step f2_std_2_3_tw$tw    PASS 1800 cert_f2.py std 2 3+ tw=$tw
  step f2_std_3_1_tw$tw    PASS 1800 cert_f2.py std 3+ 1 tw=$tw
  step f2_std_3_2_tw$tw    PASS 1800 cert_f2.py std 3+ 2 tw=$tw
  step f2_std_3_3_tw$tw    PASS 3600 cert_f2.py std 3+ 3+ tw=$tw
  step f2_two_4_1_tw$tw    PASS 1800 cert_f2.py two 4+ 1 tw=$tw
  step f2_two_4_2_tw$tw    PASS 1800 cert_f2.py two 4+ 2 tw=$tw
  step f2_two_4_3_tw$tw    PASS 1800 cert_f2.py two 4+ 3+ tw=$tw
  step f2_anti_3_1_tw$tw   PASS 1800 cert_f2.py anti 3+ 1 tw=$tw
  step f2_anti_3_2_tw$tw   PASS 1800 cert_f2.py anti 3+ 2 tw=$tw
  step f2_anti_3_3_tw$tw   PASS 1800 cert_f2.py anti 3+ 3+ tw=$tw
done
step f2_mix_2_2      PASS 1800 cert_f2.py mix 2 2 tw=0
step f2_mix_2_3      PASS 1800 cert_f2.py mix 2 3+ tw=0
step f2_mix_3_2      PASS 1800 cert_f2.py mix 3+ 2 tw=0
step f2_mix_3_3      PASS 1800 cert_f2.py mix 3+ 3+ tw=0
step f2_inv_1_3      PASS 1800 cert_f2.py inv 1 3+ tw=0
step f2_inv_3_1      PASS 1800 cert_f2.py inv 3+ 1 tw=0
step f2_inv_2_2      PASS 1800 cert_f2.py inv 2 2 tw=0
step f2_inv_2_3      PASS 1800 cert_f2.py inv 2 3+ tw=0
step f2_inv_3_2      PASS 1800 cert_f2.py inv 3+ 2 tw=0
for m in $(seq 3 "$M2"); do
  step f2_inv_${m}_3p PASS 900 cert_f2.py inv "$m" 3+ tw=0
  step f2_inv_3p_${m} PASS 900 cert_f2.py inv 3+ "$m" tw=0
done
step f2_neg_points   PASS 1800 neg_f2.py
# ---- coverage of the certified statements (Theorems F1b and F2b) by the regimes above
step coverage        PASS 600  coverage.py "$M2" M1="$M1"
t_end=$(date +%s)
echo "run_all.sh finished $(date), $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "NOT ALL PASS"; fi
