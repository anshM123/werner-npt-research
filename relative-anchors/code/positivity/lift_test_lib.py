"""Helpers shared by col_lift.py: the blocks of D' (as in ylib.Builder.Dinv, before inversion) and an exact PSD test
at a rational point."""
from ylib import Mel, Nv, halfex, wcol


def Dprime_blocks(Bd, k):
    R = Bd.R
    out = {}
    blocks = {}
    for y in Bd.cols:
        blocks.setdefault(wcol(y), []).append(y)
    for key, ys in blocks.items():
        Bm = []
        for y in ys:
            row = []
            for yp in ys:
                v = R.zero
                for t in range(k):
                    for tp in range(k):
                        mv = Mel(Bd.e[t], y, Bd.e[tp], yp)
                        if mv:
                            v += R.one * mv * R.mono(halfex(Nv(y), Nv(yp), Bd.e[t], Bd.e[tp]))
                row.append(v)
            Bm.append(row)
        out[key] = (ys, Bm)
    return out


def psd_point(M):
    n = len(M)
    A = [r[:] for r in M]
    act = list(range(n))
    while True:
        cand = [i for i in act if A[i][i] > 0]
        if not cand:
            break
        i = cand[0]
        rest = [r for r in act if r != i]
        for r in rest:
            if A[r][i] != 0:
                f = A[r][i] / A[i][i]
                for l in rest:
                    A[r][l] -= f * A[i][l]
        act = rest
    return all(A[i][i] == 0 for i in act) and all(A[i][j] == 0 for i in act for j in act)
