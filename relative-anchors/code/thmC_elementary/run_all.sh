#!/usr/bin/env bash
# Every check of NOTES.md, with logs in logs/.  Needs Python with NumPy, SymPy and python-flint.
cd "$(dirname "$0")"
mkdir -p logs
for s in symmetric_identity antisymmetric_bound margin scalars; do
  python $s.py > logs/$s.log 2>&1; echo "$s: $(tail -n 1 logs/$s.log)"
done
rm -rf __pycache__
