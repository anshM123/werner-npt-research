#!/usr/bin/env bash
# Reruns, under guard_big.sh (5.2 GB), every certificate of run_k8.sh that guard.sh stopped for memory, once at least
# 8 GB of RAM is free, until run_k8.sh reaches its coverage step.  Writes logs/fix_big.out.
cd "$(dirname "$0")"
free_gb() { powershell -NoProfile -Command "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB)" 2>/dev/null | tr -d '\r'; }
while ! grep -q "coverage_k8" logs/run_k8.out; do
  for log in $(grep -l "guard: memory" logs/run_cert_k8_*.log 2>/dev/null); do
    name=$(basename "$log" .log); name=${name#run_cert_k8_}
    [ -f "logs/rerun_done_$name" ] && continue
    while [ "$(free_gb)" -lt 8 ]; do sleep 60; done
    f=${name%%_*}; rest=${name#*_}; a=${rest%_*}; b=${rest#*_}
    echo "$(date): rerunning $f $a $b with 5.2 GB (free $(free_gb) GB)" >> logs/fix_big.out
    cp "$log" "logs/first_attempt_run_cert_k8_$name.log"
    bash guard_big.sh "$log" 7200 cert_sym.py "$f" "$a" "$b" k=8 "save=logs/cert_k8_${f}_${a}_${b}.pkl"
    echo "$(date): $name -> $(tail -n 1 "$log" | tr -d '\r')" >> logs/fix_big.out
    touch "logs/rerun_done_$name"
  done
  sleep 120
done
echo "$(date): coverage step reached, watcher stops" >> logs/fix_big.out
