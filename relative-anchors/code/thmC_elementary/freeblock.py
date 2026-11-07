"""The block of Theorem C on (X_U, q) from the closed form ../typeI_G_u3/cf.py, exactly (python-flint) or in floating
point, and the form F(X) = min_q of it on X in R^{U x U} (row-major vec(X))."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'typeI_G_u3'))
from cf import TypeI


def free_block(wT, wU, F):
    """Returns the block C on z = (X_U row-major, q) as a list of lists of numbers of type F."""
    w = [F(x) for x in list(wT) + list(wU)]
    m, u = len(wT), len(wU)
    T = list(range(m)); U = list(range(m, m + u))
    ti = TypeI(w, T, F=F)
    M = ti.build()
    zi = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    return [[M[i][j] for j in zi] for i in zi]


def q_eliminated(C):
    """F(X) = C_xx - c_xq c_qx / C_qq, as a list of lists."""
    n = len(C) - 1
    cqq = C[n][n]
    return [[C[i][j] - C[i][n] * C[n][j] / cqq for j in range(n)] for i in range(n)]
