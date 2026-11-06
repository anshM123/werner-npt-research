#!/usr/bin/env bash
# peakrun.sh LOG SECONDS args...: guard.sh with the peak memory (sampled every second) appended to LOG
log="$1"; lim="$2"; shift 2
bash ./guard.sh "$log" "$lim" "$@" &
g=$!
peak=0
sleep 1
while kill -0 $g 2>/dev/null; do
  for w in $(ps -W | awk -v gp=$g '{print $1, $2, $4}' | awk '{print $3}'); do :; done
  w=$(ps -ef | awk -v gp=$g '$3==gp {print $2}' | head -1)
  if [ -n "$w" ]; then
    wp=$(ps -W | awk -v p=$w '$1==p {print $4}')
    kb=$(tasklist //FI "PID eq $wp" //FO CSV //NH 2>/dev/null | awk -F'","' '{gsub(/[^0-9]/,"",$5); print $5}')
    [ -n "$kb" ] && [ "$kb" -gt "$peak" ] 2>/dev/null && peak=$kb
  fi
  sleep 1
done
wait $g
st=$?
echo "peak memory ${peak} KB" >> "$log.peak"
exit $st
