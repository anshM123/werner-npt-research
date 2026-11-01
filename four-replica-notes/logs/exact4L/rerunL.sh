#!/bin/bash
# full rerun of the exact4L chain in a separate folder (junctions ../exact4, ../exact4m, ../rep4 -> the real folders)
date
python peakmem.py python build_linL.py
python peakmem.py python facedataL.py | grep -v "^   entry"
python peakmem.py python spanL.py
python peakmem.py python certL.py | grep -vE "^\s+[0-9]+\s"
python peakmem.py python verifyL.py
md5sum linL.pkl facedataL.pkl facedataLr.pkl float_dualL.pkl certL.pkl
date
