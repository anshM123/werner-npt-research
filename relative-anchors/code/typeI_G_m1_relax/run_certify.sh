#!/usr/bin/env bash
# one Python process at a time, each under guard.sh (1.0 GB cap)
cd "$(dirname "$0")"
bash guard.sh logs_certify_small.log 900 certify_points.py 3,4,5,6,7,8,9,10 24 11
echo "small exit $?" > logs_certify_status.txt
bash guard.sh logs_certify_large.log 2400 certify_points.py 12,16,20,24,32 8,8,4,4,2 13 || true
echo "large exit $?" >> logs_certify_status.txt
