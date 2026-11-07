#!/usr/bin/env bash
# development helper (not part of any proof): adds workers for the cheap MU families while fewer than six Python
# processes of these runs are running (counted from the PID files of the family logs), then, when every family of
# the two MU lists has finished, runs the coverage
cd "$(dirname "$0")/.."
tagof() { echo "$1" | sed 's/-/m/g; s/[\/,]/_/g'; }
running() {
  n=0
  for f in logs/sd_*.log.pid; do
    p=$(awk '{print $1}' "$f")
    if [ -n "$p" ] && tasklist //FI "PID eq $p" //NH 2>/dev/null | grep -qi python; then n=$((n + 1)); fi
  done
  echo $n
}
unclaimed() {
  n=0
  for f in $(tr -d '\r' < lists/MU_sym_2e_cheap.txt); do [ -d "claims/sd_$(tagof "$f")" ] || n=$((n + 1)); done
  echo $n
}
finished() {
  for f in $(tr -d '\r' < lists/MU_sym_2e_cheap.txt) $(tr -d '\r' < lists/MU_sym_2d_fams.txt); do
    l="logs/sd_$(tagof "$f").log"
    [ -f "$l" ] || return 1
    last=$(tail -n 1 "$l" | tr -d '\r')
    case "$last" in PASS|FAIL|guard:*) ;; *) return 1 ;; esac
  done
  return 0
}
j=3
while [ "$(unclaimed)" -gt 0 ]; do
  if [ "$(running)" -lt 6 ]; then
    echo "$(date) $(running) running, $(unclaimed) unclaimed: one more worker (j$j)"
    JOIN=1 NW=1 bash run_sd.sh lists/MU_sym_2e_cheap.txt > "logs/run_sd_MU2e_cheap_j$j.out" 2>&1 &
    j=$((j + 1))
    sleep 30
  fi
  sleep 10
done
until finished; do sleep 15; done
echo "$(date) every MU family finished: coverage"
bash guard.sh logs/run_coverage_st4.log 3600 coverage_st3.py 120 29 30 res/*.pkl res3/q3.pkl
echo "$(date) coverage: $(tail -n 1 logs/run_coverage_st4.log)"
wait
echo "$(date) all done"
