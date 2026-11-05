"""Check: closed form (mid.Mid.full, units 2^Z = 4) against the direct construction, k = 4, 5 (float),
and the claimed kernel (p (x) H_{k-2} and vac (x) J^-|1..1>) exactly at rational weights, k = 4, 5, 6."""
import numpy as np
from fractions import Fraction as Fr
from mid import Mid
from direct import direct_block

ok = True
rng = np.random.default_rng(1)
for k in (4, 5):
    for trial in range(3):
        w = list(rng.random(k) + 0.2)
        if trial == 2:
            w[0] = 0.05
        Yd = direct_block(w)
        Yc = np.array(Mid(w, float).full(), float) / 4
        err = np.abs(Yd - Yc).max()
        ev = np.linalg.eigvalsh(Yd)
        print(f"k={k} w={np.round(w, 3)}: |direct - closed| = {err:.2e}; eig min {ev[0]:+.2e}, #|ev|<1e-9 = {np.sum(np.abs(ev) < 1e-9)}")
        ok &= err < 1e-11

for k in (4, 5, 6):
    w = [Fr(int(x), 7) for x in rng.integers(2, 25, size=k)]
    m = Mid(w)
    Y = m.full()
    st = m.states()
    N = len(st)
    vecs = []
    for P in m.P:
        v = [Fr(0)] * N
        for i, X in enumerate(st):
            if X[0] == 'P' and X[2] == P:
                v[i] = m.w[X[1]]
        vecs.append(v)
    v = [Fr(0)] * N
    for i, X in enumerate(st):
        if X[0] == 'Q':
            v[i] = m.w[X[1]]
    vecs.append(v)
    res = max(abs(sum(Y[i][j] * vv[j] for j in range(N))) for vv in vecs for i in range(N))
    print(f"k={k} exact: |Y kappa| over the {len(vecs)} kernel vectors = {res}")
    ok &= res == 0
print("PASS" if ok else "FAIL")
