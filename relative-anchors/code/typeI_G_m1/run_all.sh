#!/usr/bin/env bash
# Reruns every step of this folder, one Python process at a time, each under guard.sh (killed above 2.8 GB or after
# its time limit).  Prints PASS/FAIL per step and ALL PASS at the end.  Logs: logs/run_<step>.log, summary in
# logs/run_all.out when called as in NOTES.md.  About 20 minutes.  (The symbolic three-level reduction
# three_level_fast.py does not finish in reasonable time and is not run here; NOTES.md, Section 6.)
cd "$(dirname "$0")"
mkdir -p logs
fail=0
step() {
  name="$1"; lim="$2"; shift 2
  ./guard.sh "logs/run_${name}.log" "$lim" "$@"
  st=$?
  if [ $st -eq 0 ] && tail -1 "logs/run_${name}.log" | grep -q '^PASS'; then
    echo "PASS  $name"
  else
    echo "FAIL  $name (exit $st; see logs/run_${name}.log)"
    fail=1
  fi
}
# 1. exact closed forms for u = 1, 2, 3 in the gauge x_t = 0; Conjecture M there by coefficient signs
step closed_t_u1 120 closed_t.py 1
step closed_t_u2 120 closed_t.py 2
step closed_t_u3 600 closed_t.py 3
step check_m_small 300 check_m_small.py
# 2. equal U-weights, every u (symbolic in u), and the limit u -> infinity
step equal_u 900 equal_u.py
step check_equal_u 900 check_equal_u.py
step equal_u_limit 600 equal_u_limit.py
# 3. two-level U-weights, every (u1, u2) (symbolic in u1, u2)
step two_level_inv_one 600 two_level_fast.py inv one
step two_level_inv_gen 1200 two_level_fast.py inv gen
step two_level_std_u1is2 600 two_level_fast.py std u1is2
step two_level_std_gen 600 two_level_fast.py std gen
step check_two_level_one 600 check_two_level.py one
step check_two_level_gen 1200 check_two_level.py gen
step check_two_level_std 1200 check_two_level_std.py
# 3b. second, independent implementation of the invariant two-level reduction (case one) and exact comparison
step two_level_one_exprtree 600 two_level.py one
step cmp_impl 120 cmp_impl.py
# 4. exact checks: rational points for u <= 10; three levels at rational weights, uniformly in u3; negative controls
step check_points 2400 check_points.py 10 28
step three_level_points 1200 three_level_points.py 12
step where_min 900 where_min.py 8
step controls 600 controls.py
# 5. floating-point exploration (no claim rests on it)
step ex_sep_run 300 ex_sep_run.py
step ex_ev 600 ex_ev.py
if [ $fail -eq 0 ]; then echo "ALL PASS"; else echo "SOME STEPS FAILED"; fi
