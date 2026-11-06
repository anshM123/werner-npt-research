#!/usr/bin/env bash
# Reruns every check and certificate of this folder, one Python process at a time, each under guard.sh (2.5 GB, time
# limit per step).  PASS or FAIL per step, ALL PASS at the end.  Steps marked "expect FAIL" are negative controls of
# certificates: they pass if the certificate run ends with FAIL (and the process was not stopped by the guard).
# Usage: bash run_all.sh > logs/run_all.out   (KMAX = largest k of the fixed-k certificates, default 32;
#                                              W = width of the bands of j at both ends, default 8)
cd "$(dirname "$0")"
mkdir -p logs
all=1
t_start=$(date +%s)
KMAX=${KMAX:-32}
W=${W:-8}

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

echo "run_all.sh started $(date), KMAX = $KMAX, W = $W"
# ---- the reduction: orbit calculus against the Lemma B builder; symbolic against numeric; window against all nodes
step builder          PASS 1800 check_builder.py 6 6 16
# ---- completeness of the reduction (exact dimension counts; spectra in floating point)
step spectrum         PASS 1800 check_spectrum.py 6
# ---- the fixed-k certificates, k = 2..KMAX, every j, every i (saves the certified forms)
step cert_fixed       PASS 7200 cert_fixed.py 2 "$KMAX" 1 logs/cert_fixed.pkl
# ---- coverage of every (k, j, i), and the saved certified forms against the Lemma B builder
step coverage         PASS 1800 coverage.py logs/cert_fixed.pkl 2 "$KMAX" 6
# ---- negative controls: B -> (1 + 1/1000)B
step neg_cert         FAIL 1800 cert_fixed.py 2 16 1001/1000
step neg_points       PASS 1800 neg_points.py 6
# ---- families symbolic in the sizes, sigma free: every family with one symbolic size; the bands of width W at both
#      ends of j (symbolic in k); their coverage; negative control
step fam1             PASS 3600 prove_fam.py 1 1
step band_bottom      PASS 7200 prove_band.py bottom "$W"
step band_top         PASS 7200 prove_band.py top "$W"
step coverage_fam     PASS 1800 coverage_fam.py 60 "$KMAX" "band:$W:$W" logs/fam_1_1.pkl "logs/band_bottom_$W.pkl" "logs/band_top_$W.pkl"
step neg_fam          FAIL 3600 prove_fam.py 1 1 1001/1000
# ---- the invariant part (i = 0) for every k and every j: the family [q0 >= 4, t >= 4, i = 0] and its coverage
step fam_invariant    PASS 1800 prove_lowirr.py 0
step coverage_inv     PASS 1800 coverage_irr.py 60 "$KMAX" 0 logs/fam_lowirr_0.pkl logs/fam_1_1.pkl "logs/band_bottom_$W.pkl" "logs/band_top_$W.pkl"
t_end=$(date +%s)
echo "total $((t_end - t_start)) s"
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
