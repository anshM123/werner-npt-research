#!/usr/bin/env bash
# guard.sh LOG SECONDS args...   runs one Python process (python args...), writes its output to LOG, and kills it by
# its own Windows PID if its memory exceeds 2.5 GB or if it runs longer than SECONDS (a line "guard: ..." is then
# appended to LOG and the exit status is 3).  The Windows PID and the arguments are written to LOG.pid.
log="$1"; shift
lim="$1"; shift
python "$@" > "$log" 2>&1 &
pid=$!
sleep 0.5
wpid=$(ps -p "$pid" 2>/dev/null | awk 'NR==2 {print $4}')
echo "$wpid $*" > "$log.pid"
start=$(date +%s)
killed=0
while kill -0 "$pid" 2>/dev/null; do
  if [ -n "$wpid" ]; then
    kb=$(tasklist //FI "PID eq $wpid" //FO CSV //NH 2>/dev/null | awk -F'","' '{gsub(/[^0-9]/,"",$5); print $5}')
    if [ -n "$kb" ] && [ "$kb" -gt 2500000 ] 2>/dev/null; then
      echo "guard: memory ${kb} KB above 2.5 GB, killed PID $wpid" >> "$log"
      taskkill //PID "$wpid" //F > /dev/null 2>&1
      killed=1
      break
    fi
    now=$(date +%s)
    if [ $((now - start)) -gt "$lim" ]; then
      echo "guard: time limit ${lim} s exceeded, killed PID $wpid" >> "$log"
      taskkill //PID "$wpid" //F > /dev/null 2>&1
      killed=1
      break
    fi
  fi
  sleep 1
done
wait "$pid" 2>/dev/null
st=$?
if [ "$killed" = 1 ]; then exit 3; fi
exit $st
