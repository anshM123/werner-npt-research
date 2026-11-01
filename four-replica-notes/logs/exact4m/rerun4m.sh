#!/bin/bash
# full rerun of the chain after the points (reproducibility check); run from exact4m
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
date
python peakmem.py python build_lin4m.py | grep -v "identity ORTH\|identity MARG"
python peakmem.py python face4m.py z_k2_a z_k2_b z_k4_a z_k4_b z_in_a z_in_b z_in_c z_in_d z_sf_a z_sf_b z_sf_c z_sf_d z_sf_e z_sf_f z_sf_g z_sf_h | grep -v "^   entry"
python peakmem.py python facedata4m.py | grep -v "^   entry"
python peakmem.py python span4m.py
python peakmem.py python cert4m.py
python peakmem.py python verify4m.py p_ent p_ent2
md5sum lin4m.pkl face4m.pkl facedata4m.pkl facedata4mr.pkl float_dual4m.pkl cert4m.pkl
date
