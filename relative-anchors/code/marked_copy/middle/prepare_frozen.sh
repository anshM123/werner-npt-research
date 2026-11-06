#!/usr/bin/env bash
# The runs here read the outcomes of the band families and of Theorems F and I of the parent folder from logs/frozen/.
# Run the parent's run_all.sh first (it writes ../logs/*.pkl), then this script, then run_all.sh and the extensions.
cd "$(dirname "$0")"
mkdir -p logs/frozen
for f in band_bottom_8.pkl band_top_8.pkl fam_1_1.pkl fam_lowirr_0.pkl run_all.out; do cp "../logs/$f" logs/frozen/; done
echo "copied the outcomes of ../run_all.sh into logs/frozen/"
