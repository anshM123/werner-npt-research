#!/usr/bin/env bash
# Every step of NOTES.md, logs in logs/.  Needs Node (no packages) and Python with python-flint; the Python checks use
# ../typeI_G_m1_relax/psi_exact.py.  One process at a time; the certificate takes a few minutes.
cd "$(dirname "$0")"
mkdir -p logs
all=1
run() { local name="$1"; shift; "$@" > "logs/$name.log" 2>&1; local last=$(tail -n 1 "logs/$name.log" | tr -d '\r'); echo "$name: $last"; }
python run_js.py one_letter_certificate.js > logs/one_letter_certificate.log 2>&1
grep -q '"negativeCoefficients": 0' logs/one_letter_certificate.log && grep -q '"positiveCoefficients": 20406' logs/one_letter_certificate.log && echo "certificate: PASS" || { echo "certificate: FAIL"; all=0; }
python run_js.py relax_comparison.js sfQReviewA510RelaxAllVerify > logs/relax_comparison.log 2>&1
[ "$(grep -c '"entrywiseIdentity": true' logs/relax_comparison.log)" = 12 ] && echo "relax_comparison: PASS" || { echo "relax_comparison: FAIL"; all=0; }
run stress_profiles python stress_profiles.py
run crossbase_s5 python crossbase.py 5 16
run bridge_s1 python bridge.py 1 6
run test_A30_s1 python test_A30.py 1 12
run test_A30_control_half python test_A30.py 1 8 1 2     # expected FAIL (negative control)
rm -f cross_cases.json cross_harness.js
for f in stress_profiles crossbase_s5 bridge_s1 test_A30_s1; do [ "$(tail -n 1 logs/$f.log | tr -d '\r')" = PASS ] || all=0; done
[ "$(tail -n 1 logs/test_A30_control_half.log | tr -d '\r')" = FAIL ] || all=0
if [ "$all" = 1 ]; then echo "ALL PASS"; else echo "SOME STEP FAILED"; fi
