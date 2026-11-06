#!/usr/bin/env bash
# The regime of Theorem F2b left open by run_all.sh and run_fext.sh: family 2 (T = {t}), invariant part, both U-levels
# symbolic (3+, 3+), t-weight 0, by determinant and inertia (f2det.py, Section 5d of NOTES.md): bounds and candidate
# denominator, the slices x1 = 0..B (NW workers in parallel, default 1; the logged run used 3), the interpolation in x1
# with the signs, the base point; then the coverage of families 1 and 2 with both open regimes included.
# Each Python process runs under guard.sh (2.5 GB).  Data in f2det_out/ (not kept in the repository).
# Usage: NW=3 bash run_open2.sh > logs/run_open2.out
cd "$(dirname "$0")"
mkdir -p logs
NW=${NW:-1}
allok=1
step() {
  name="$1"; exp="$2"; lim="$3"; shift 3
  t0=$(date +%s)
  bash guard.sh "logs/run_${name}.log" "$lim" "$@"
  last=$(tail -n 1 "logs/run_${name}.log" | tr -d '\r')
  if [ "$last" = "$exp" ]; then r=PASS; else r=FAIL; allok=0; fi
  echo "$r  $name  (python $*; expect $exp; $(( $(date +%s) - t0 )) s)"
}
echo "run_open2.sh started $(date)"
step f2det_bounds PASS 3600 f2det.py bounds
t0=$(date +%s)
for w in $(seq 0 $((NW - 1))); do
  bash guard.sh logs/run_f2det_slices_$w.log 108000 f2det.py slices $w $NW &
done
wait
for w in $(seq 0 $((NW - 1))); do
  last=$(tail -n 1 logs/run_f2det_slices_$w.log | tr -d '\r')
  if [ "$last" = PASS ]; then r=PASS; else r=FAIL; allok=0; fi
  echo "$r  f2det_slices_$w  (python f2det.py slices $w $NW; expect PASS; $(( $(date +%s) - t0 )) s)"
done
step f2det_interp PASS 7200 f2det.py interp
step f2det_final PASS 1800 f2det.py final
step coverage_open2 PASS 1800 coverage.py 40 M1=40 KMAX=120 F1OPEN=1 F2OPEN=1
# checks, not needed by the proof: the interpolant at points outside the grid; the negative control (expected FAIL)
step f2check_points PASS 3600 f2check.py points
step f2check_neg FAIL 3600 f2check.py neg
if [ "$allok" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
