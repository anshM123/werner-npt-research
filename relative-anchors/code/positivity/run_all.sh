#!/usr/bin/env bash
# Reruns every step of this folder, one Python process at a time, each under guard.sh (killed above 2.8 GB).
# Prints PASS or FAIL per step and ALL PASS at the end.  A step passes if its log ends with PASS and contains the
# expected result line.  LONG=1 adds the 1023-minor run (about 7 min) and the random k = 4 sample (about 5 min).
cd "$(dirname "$0")"
mkdir -p logs
ALL=0
step () {
  name="$1"; expect="$2"; shift 2
  log="logs/run_${name}.log"
  ./guard.sh "$log" "$@"
  if grep -q "^PASS" "$log" && ! grep -q "^guard:" "$log" && { [ -z "$expect" ] || grep -qF "$expect" "$log"; }; then
    echo "PASS  $name"
  else
    echo "FAIL  $name  (see $log)"; ALL=1
  fi
}
if [ ! -f blocks_sym_k4.pkl ]; then
  (cd ../middle_k4 && python build_sym.py > /dev/null 2>&1) && cp ../middle_k4/blocks_sym.pkl blocks_sym_k4.pkl
fi
step build_k2 "" build_blocks.py 2
step build_k3 "" build_blocks.py 3
step allminors_k2 "TOTAL: 46 principal minors, positive cores 46, mixed 0, point-check failures 0" allminors.py 2
step allminors_k3 "TOTAL: 916 principal minors, positive cores 916, mixed 0, point-check failures 0" allminors.py 3
step allminors_k4_small "TOTAL: 22 principal minors, positive cores 22, mixed 0, point-check failures 0" allminors.py 4 all 4
step rowgroup_minors_k4 "TOTAL: 447 principal minors of row groups, positive cores 447, mixed 0, point-check failures 0" rowgroup_minors.py 4
step mono_psd_k2 "RESULT k = 2: 10 of 12 blocks" mono_psd.py 2
step mono_psd_k3 "RESULT k = 3: 13 of 20 blocks" mono_psd.py 3
step mono_psd_k4_small "RESULT k = 4: 2 of 2 blocks" mono_psd.py 4 4
step newton_k2 "RESULT k = 2: 1 of 12 blocks have an exact obstruction, 1 of them also" newton_obstr.py 2
step newton_k3 "RESULT k = 3: 6 of 20 blocks have an exact obstruction, 2 of them also" newton_obstr.py 3
step newton_k4 "RESULT k = 4: 5 of 7 blocks have an exact obstruction, 3 of them also" newton_obstr.py 4 nolp
step obstr_rows_k3 "RESULT k = 3: inside row groups 26, between row groups 27" obstr_rows.py 3
step obstr_rows_k4 "RESULT k = 4: inside row groups 124, between row groups 200" obstr_rows.py 4
step obstr_cols_k3 "RESULT k = 3: same row 26, same column 0, neither 27" obstr_cols.py 3
step obstr_cols_k4 "RESULT k = 4: same row 124, same column 0, neither 200" obstr_cols.py 4
step col_lift_k3 "RESULT k = 3: 14 column groups with at least 2 states; lift exact and PSD for 14; a fixed-vector representation found for 12" col_lift.py 3
step col_lift_k4 "RESULT k = 4: 15 column groups with at least 2 states; lift exact and PSD for 15; a fixed-vector representation found for 13" col_lift.py 4
step farkas_k2 "RESULT k = 2: 1 of 1 obstructed blocks have an exact certificate" rescale_farkas.py 2
step farkas_k3 "RESULT k = 3: 5 of 6 obstructed blocks have an exact certificate" rescale_farkas.py 3
step farkas_k4 "RESULT k = 4: 5 of 5 obstructed blocks have an exact certificate" rescale_farkas.py 4
step rescale_k2 "" mono_rescale.py 2
step rescale_k3 "" mono_rescale.py 3
step inv_k2 "" inv_test.py 2
step inv_k3 "" inv_test.py 3
step jet_k2 "RESULT k = 2: 0 of 2 blocks" jet_obstr.py 2
step jet_k3 "RESULT k = 3: 0 of 7 blocks" jet_obstr.py 3
step polya_k2 "RESULT k2_j1_(1, 1, 0): no representation found" polya_mult.py 2 12 "k2_j1_(1, 1, 0)"
step polya_k3 "RESULT k3_j2_(1, 1, 1, 0): no representation found" polya_mult.py 3 6 "k3_j2_(1, 1, 1, 0)"
step linear_mult "best (3, 1, 1, 1)" linear_mult.py
step k2_rescaled "" k2_rescaled.py
step rowgroup_identity "" rowgroup_identity.py
step rowgroup_lift "k = 4: 21368 entries checked so far, all equal: True" rowgroup_lift.py
step lift_k2_j0 "after the LP rescaling: 0 obstructed pairs; monomial-PSD: 1 non-PSD of 7" lift_test.py 2 "k2_j0_(1, 0, 0)"
step lift_k3_j0_L1 "monomial rescalings of the lift: certificate -8/13" lift_test.py 3 "k3_j0_(1, 0, 0, 0)"
step lift_k3_j1_I "monomial rescalings of the lift: certificate -2/3" lift_test.py 3 "k3_j1_(1, 1, 0, 0)"
step lift_k3_j1_L1 "monomial rescalings of the lift: certificate -40/69" lift_test.py 3 "k3_j1_(2, 0, 0, 0)"
step lift_k3_j1_T01 "after the LP rescaling: 0 obstructed pairs; monomial-PSD: 12 non-PSD of 44" lift_test.py 3 "k3_j1_(0, 1, 1, 0)"
if [ "${LONG:-0}" = "1" ]; then
  step allminors_k4_j1_T0 "positive cores 1023, mixed 0, point-check failures 0" allminors.py 4 k4_j1_T0
  step randsub_k4 "TOTAL: 504 principal minors, positive cores 504, mixed 0, point-check failures 0" allminors.py 4 all 99 sub:12:7:2026
fi
if [ $ALL = 0 ]; then echo "ALL PASS"; else echo "SOME STEPS FAILED"; fi
