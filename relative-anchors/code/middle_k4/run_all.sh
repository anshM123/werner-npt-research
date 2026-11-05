#!/bin/sh
# (R_inf)_1 and (R_inf)_2 at k = 4: letter-1 and type I torus blocks.  Reruns every step, one Python process at a
# time, prints PASS/FAIL per step and ALL PASS / SOME FAIL at the end.  Logs in logs/.  Needs Python 3 with numpy and
# python-flint (0.9 used).  ylib.py and flintfield.py are copies from relative-anchors/code/.
# LONG=0 skips the three large certificates (sizes 18, 21, 18; about 10, 11 and 65 minutes); see NOTES.md.
cd "$(dirname "$0")"
mkdir -p logs
: > logs/run_all.log
fail=0
run() {
  name=$1; shift
  tag=$(echo "$name $*" | sed 's/ *$//' | tr ' ' '_')
  start=$(date +%s)
  python -u ${name}.py "$@" > logs/${tag}.log 2>&1
  last=$(tail -n 1 logs/${tag}.log | tr -d '\r')
  end=$(date +%s)
  echo "$name $*: $last ($((end - start)) s)" | tee -a logs/run_all.log
  [ "$last" = "PASS" ] || fail=1
}
run num_k4 30            # torus structure (exact), known kernels and gauge-fixed positivity at 33 rational points (exact LDL)
run build_sym            # exact blocks over Q(s1, s2, s3), s0 = 1 -> blocks_sym.pkl; known kernels annihilated identically
run check_direct         # blocks = dense construction from the definitions (float); symbolic = exact point builder (exact)
# ---- certificates: all leading principal minors of the gauge-fixed block have only positive coefficients
run bareiss_cert 1 01 nat     # j = 1, T = {t, t'}   (size 3)
run bareiss_cert 2 012 nat    # j = 2, T = {t, t', t''} (size 4)
run bareiss_cert 1 0 nat      # j = 1, T = {t}        (size 10, about 20 s)
run bareiss_cert 2 01 nat     # j = 2, T = {t, t'}    (size 11, about 2.5 min)
run elim_minors 1 01 nat      # independent route (elimination over the rational function field), small blocks
run elim_minors 2 012 nat
if [ "${LONG:-1}" = "1" ]; then
  run bareiss_cert 1 - nat    # j = 1, letter-1 block (size 18, about 10 min)
  run bareiss_cert 2 - nat    # j = 2, letter-1 block (size 21); see NOTES.md for the run time
  run bareiss_cert 2 0 nat    # j = 2, T = {t}        (size 18, about 65 min)
fi
# end-to-end: certified minors = leading minors of the directly built block at rational points (uses the cert_*.pkl
# files written above; with LONG=0 the large ones come from earlier runs of the same commands)
run check_minors_points   # all six certificates (default argument list in the script)
run controls             # negative controls: (1 + 1/1000) B fails on every block, exact witnesses
if [ $fail -eq 0 ]; then echo "ALL PASS" | tee -a logs/run_all.log; else echo "SOME FAIL" | tee -a logs/run_all.log; fi
