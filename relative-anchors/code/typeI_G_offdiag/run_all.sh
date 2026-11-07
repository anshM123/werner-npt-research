#!/usr/bin/env bash
# Every check of NOTES.md, with logs in logs/.  Needs Python with python-flint, NumPy and SymPy.
cd "$(dirname "$0")"
mkdir -p logs
run() { local log="logs/$1.log"; shift; python "$@" > "$log" 2>&1; echo "$(basename "$log" .log): $(tail -n 1 "$log")"; }
run lemmaq_exact      lemmaq_exact.py 1 8
run lemmaL_check      lemmaL_check.py
run check_exact       check_exact.py 1 8
run exact_gg2_k38     exact_gg2.py 11 12 3/8
run exact_gg2_k14     exact_gg2.py 11 12 1/4
run exact_formal_k38  exact_formal.py 5 20 3/8
run exact_formal_k14  exact_formal.py 5 20 1/4
rm -rf __pycache__
