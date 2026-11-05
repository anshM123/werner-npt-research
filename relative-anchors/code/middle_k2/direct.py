"""Direct construction of the letter-1 block of (R_inf)_{k-2} from M, V_p, D (check_bath_reduction.group_ops)."""
import sys, itertools
import numpy as np
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'single-excitation-range', 'code'))
from check_bath_reduction import group_ops
from mid import pairs


def hole_string(k, holes):
    return tuple(0 if t in holes else 1 for t in range(k))


def direct_block(w, j=None):
    """Unrescaled Y on the letter-1 block of (R_inf)_{k-2}, states in the order of Mid.states()."""
    k = len(w)
    pG = {}
    for t in range(k):
        x = [0] * k
        x[t] = 1
        pG[tuple(x)] = float(w[t])
    _, Aent, Bent = group_ops(pG, k)
    st = []
    for r in range(k):
        row = tuple(1 if t == r else 0 for t in range(k))
        for P in pairs(k):
            st.append((row, hole_string(k, P)))
    vac = (0,) * k
    for h in range(k):
        st.append((vac, hole_string(k, (h,))))
    Y = np.array([[Aent(a, yp, b, y) - Bent(a, yp, b, y) for (b, y) in st] for (a, yp) in st])
    return (Y + Y.T) / 2
