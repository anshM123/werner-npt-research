#!/usr/bin/env bash
# Every check of this folder, with its log in logs/.  Needs Python with NumPy.
cd "$(dirname "$0")"
mkdir -p ../logs
python route_checks.py > ../logs/route_checks.log 2>&1; tail -n 1 ../logs/route_checks.log
python norm_bounds_check.py > ../logs/norm_bounds_check.log 2>&1; tail -n 1 ../logs/norm_bounds_check.log
python reference_obstructions_check.py > ../logs/reference_obstructions_check.log 2>&1; tail -n 1 ../logs/reference_obstructions_check.log
