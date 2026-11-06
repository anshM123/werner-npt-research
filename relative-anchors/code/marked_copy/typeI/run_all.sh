#!/usr/bin/env bash
# Every step of this folder (NOTES.md), one Python process at a time, each under guard.sh (2.5 GB, a time limit per
# step).  PASS or FAIL per step, ALL PASS at the end.  Steps marked "expect FAIL" are negative controls
# (B -> (1 + 1/1000) B): they pass if the run ends with FAIL (and the process was not stopped by the guard).
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
# ---- the reduction against the Lemma B block (FastJ1), dimension counts, kernel vectors: k <= 8, every j >= 1
step explicit_3_8     PASS 1800 check_explicit.py 3 8 1 2,1/3
# ---- symbolic forms against numeric ones, the Newton window against all nodes, the sector denominators: k <= 10
step symnum_3_10      PASS 1800 check_symnum.py 3 10
# ---- the fixed-k certificates, 3 <= j <= k - 1, every tau, both cases, both T'-modes, every i
step cert_08_12       PASS 1800 cert_t.py 8 12 1 logs/cert_08_12.pkl
step cert_13_15       PASS 1800 cert_t.py 13 15 1 logs/cert_13_15.pkl
step cert_16_17       PASS 1800 cert_t.py 16 17 1 logs/cert_16_17.pkl
step cert_18          PASS 1800 cert_t.py 18 18 1 logs/cert_18.pkl
step cert_19          PASS 1800 cert_t.py 19 19 1 logs/cert_19.pkl
step cert_20          PASS 1800 cert_t.py 20 20 1 logs/cert_20.pkl
# ---- coverage of every type I block with a letter 1, 8 <= k <= 20
step coverage_8_20    PASS 1800 coverage_t.py 8 20 logs/cert_08_12.pkl logs/cert_13_15.pkl logs/cert_16_17.pkl logs/cert_18.pkl logs/cert_19.pkl logs/cert_20.pkl
# ---- negative control of the certificates: B -> (1 + 1/1000) B (expect FAIL: some form fails at every k)
step neg_cert_8_12    FAIL 1800 cert_t.py 8 12 1001/1000 -
# ---- whole blocks end to end at points (exact, no reduction), with the scaled control
step e2e_8            PASS 1800 e2e_t.py 8 2,1/3
step e2e_9            PASS 3000 e2e_t.py 9 3/2
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
