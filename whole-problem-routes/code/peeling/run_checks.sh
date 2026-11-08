#!/bin/sh
# Reproduces the logs of code/peeling (Proposition 7.30, Remark 7.31).  A few minutes in total, small memory.
cd "$(dirname "$0")"
python product_probe.py > logs/product_probe.log 2>&1
for mode in tc peel1 peelbest vtc1 vtcbest; do
    python peel_search.py 2 9 $mode 12 1 > logs/peel_search_n2_$mode.log 2>&1
done
python peel_search.py 3 27 tc 12 4 > logs/peel_search_n3_tc.log 2>&1
python peel_search.py 3 27 vtcbest 12 7 > logs/peel_search_n3_vtcbest.log 2>&1
python peel_search.py 4 9 tc 8 11 > logs/peel_search_n4_tc.log 2>&1
{
    python vertex_level2.py 1 3 4 rebit
    python vertex_level2.py 1 3 4 qutrit
    python vertex_level2.py 2 3 4 rebit
    python vertex_level2.py 2 3 4 qutrit
} > logs/vertex_level2.log 2>&1
