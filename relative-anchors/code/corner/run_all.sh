#!/bin/sh
# Proposition C (README, Section 4): the last family of Theorem E3, by an exact determinant certificate.
# One Python process at a time; about 15 minutes (the 15 node determinants dominate; compute_nodes.py skips nodes
# already saved in nodes/); peak memory below 0.2 GB. nodes/ and core.pkl are regenerated and not kept in git.
set -e
cd "$(dirname "$0")"
mkdir -p logs nodes
python -u compute_nodes.py 2 16                > logs/compute_nodes.log 2>&1
python -u interpolate.py                       > logs/interpolate.log 2>&1
# the corner {q0 >= 2, t >= 6, i >= 2, tau >= 2}
python -u shift_check.py core.pkl 2 6 2 2      > logs/shift_check_corner.log 2>&1
python -u check_structure.py 2 6 2 2           > logs/check_structure_corner.log 2>&1
python -u check_det_points.py 2 6 2 2          > logs/check_det_points_corner.log 2>&1
# the whole generic family {q0 >= 2, t >= 4, i >= 2, tau >= 1}
python -u shift_check.py core.pkl 2 4 2 1      > logs/shift_check_generic.log 2>&1
python -u check_structure.py 2 4 2 1           > logs/check_structure_generic.log 2>&1
python -u check_det_points.py 2 4 2 1          > logs/check_det_points_generic.log 2>&1
# the irrep split: P invertible, Q = Q8 (+) Q2, kappa in ker Q8, Q2 positive definite
python -u check_split.py                       > logs/check_split.log 2>&1
grep -H "PASS\|FAIL\|NOT " logs/*.log
