#!/usr/bin/env bash
# Propositions D and D1 with the split (tsdet.py, NOTES.md Section 8), one family per Python process (guard.sh,
# 1.2 GB), NW workers that claim the families of a list in turn (atomic mkdir of claims/NAME).  Each family writes
# res/sd_NAME.pkl and logs/sd_NAME.log; a family whose log ends with PASS is skipped.
# LIST lines: 'MU/sym/a,c,n,n' (i and tau' symbolic: Proposition D) or any two-symbolic family (Proposition D1).
# Usage: NW=3 bash run_sd.sh LIST > logs/run_sd_NAME.out      (resumable: rerun; stale claims are removed at start
#        unless JOIN=1, which only adds workers to a run in progress)
cd "$(dirname "$0")"
NW=${NW:-3}
mkdir -p logs res claims
fams=($(tr -d '\r' < "$1"))
tagof() { echo "$1" | sed 's/-/m/g; s/[\/,]/_/g'; }
if [ "${JOIN:-0}" != 1 ]; then
  for f in "${fams[@]}"; do rmdir "claims/sd_$(tagof "$f")" 2>/dev/null; done
fi
worker() {
  for f in "${fams[@]}"; do
    nm=$(tagof "$f")
    log="logs/sd_$nm.log"
    if [ -f "$log" ] && [ "$(tail -n 1 "$log" | tr -d '\r')" = PASS ]; then continue; fi
    mkdir "claims/sd_$nm" 2>/dev/null || continue
    if [[ $f =~ ^(MT|MU)/sym/(-?[0-9]+,-?[0-9]+),n,n$ ]]; then
      bash guard.sh "$log" 36000 tsdet.py D "res/sd_$nm.pkl" "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
    else
      bash guard.sh "$log" 36000 tsdet.py D1 "res/sd_$nm.pkl" "$f"
    fi
    echo "$f: $(tail -n 2 "$log" | head -n 1 | tr -d '\r') / $(tail -n 1 "$log" | tr -d '\r')"
  done
}
for w in $(seq 1 "$NW"); do worker & done
wait
echo "run_sd.sh $1 done $(date)"
