#!/usr/bin/env bash
# Proposition D2 (tdet2.py, NOTES.md Section 8) for the families of a list (sym, i and tau' symbolic, exactly one of
# q0, t symbolic): prep one family at a time, then NW slice workers (each one Python process under guard.sh, 1.2 GB),
# which claim free columns of each family in turn and skip saved slices, then the finals one at a time.
# Usage: NW=3 bash run_d2.sh LIST > logs/run_d2_NAME.out      (resumable: rerun the same command)
#        JOIN=1 NW=1 bash run_d2.sh LIST                        (one more worker for a run in progress: no prep, no
#                                                               removal of claims, no finals)
cd "$(dirname "$0")"
NW=${NW:-3}
mkdir -p logs res d2
fams=($(tr -d '\r' < "$1"))
tagof() { echo "$1" | sed 's/-/m/g; s/[\/,]/_/g'; }
if [ "${JOIN:-0}" != 1 ]; then
  for f in "${fams[@]}"; do
    d=d2/$(tagof "$f")
    if [ ! -f "$d/meta.pkl" ]; then
      bash guard.sh "logs/d2_prep_$(tagof "$f").log" 3600 tdet2.py prep "$f" "$d"
      echo "prep $f: $(tail -n 1 "logs/d2_prep_$(tagof "$f").log")"
    fi
    rm -f "$d"/claim_*
  done
fi
for w in $(seq 1 "$NW"); do
  wid=$((w + ${WOFF:-0}))
  ( for f in "${fams[@]}"; do
      d=d2/$(tagof "$f")
      [ -f "$d/meta.pkl" ] && bash guard.sh "logs/d2_w${wid}_$(tagof "$f").log" 36000 tdet2.py slices "$d" "$wid" 0
    done ) &
done
wait
if [ "${JOIN:-0}" != 1 ]; then
  for f in "${fams[@]}"; do
    d=d2/$(tagof "$f")
    [ -f "$d/meta.pkl" ] || continue
    bash guard.sh "logs/d2_final_$(tagof "$f").log" 7200 tdet2.py final "$d" "res/det2_$(tagof "$f").pkl"
    echo "final $f: $(tail -n 1 "logs/d2_final_$(tagof "$f").log")"
  done
fi
echo "run_d2.sh $1 done $(date)"
