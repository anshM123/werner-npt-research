#!/bin/sh
# Proposition U1 (README, Section 3): type I blocks of (R_inf)_{k-1} with |U| = 1, every k, all weights.
# One Python process at a time (each well under 0.5 GB); about 25 minutes, most of it in the two symbolic identity
# checks.  gstruct.pkl, polys.pkl and polys.txt are regenerated and not kept in git.
cd "$(dirname "$0")"
mkdir -p logs
: > logs/run_all.log
run() {
  name=$1; log=$2; shift 2
  start=$(date +%s)
  python -u ${name}.py "$@" > logs/${log}.log 2>&1
  end=$(date +%s)
  echo "$name $*: $(tail -n 1 logs/${log}.log) ($((end - start)) s)" | tee -a logs/run_all.log
}
run gstruct gstruct
run test_bern test_bern
run check_closed_forms check_closed_forms
run sym_check_field sym_check_field 1 2 3
run polys polys
run check_endtoend check_endtoend
run prove_AB prove_AB
run controls controls
run sanity_float sanity_float
run sym_check sym_check 3
run sym_check_field sym_check_field_m4 4
