#!/usr/bin/env bash
# development driver: direct certificates at k = 6 for the representative configurations, smallest families first
cd "$(dirname "$0")"
for f in "$@"; do
  for c in $(python -c "import k5lib as K; j,m=K.fam('$f'); print(' '.join(f'{a},{b}' for a,b in K.rep_configs(6,m)))"); do
    a=${c%,*}; b=${c#*,}
    t0=$(date +%s)
    bash guard.sh logs/dev/cert_k6_${f}_${a}_${b}.log 2400 cert5.py $f $a $b k=6 save=logs/dev/cert_k6_${f}_${a}_${b}.pkl
    st=$?
    t1=$(date +%s)
    echo "$f $a $b: exit $st, $(tail -n1 logs/dev/cert_k6_${f}_${a}_${b}.log | tr -d '\r'), $(grep -o 'nonnegative coefficients: [A-Za-z]*; largest degrees.*' logs/dev/cert_k6_${f}_${a}_${b}.log), $((t1 - t0)) s"
  done
done
echo DRIVER DONE
