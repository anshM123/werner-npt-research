import numpy as np
from numlib import block
from forms import Forms
rng = np.random.default_rng(3)
mx = 0
for it in range(20):
    u = int(rng.integers(2, 7))
    st = float(np.exp(rng.uniform(-2, 2))); sU = np.exp(rng.uniform(-2, 2, u))
    P, lab, ne = block(st, sU, scaled=False)
    nX = u * u; iq = nX
    keep = [i for i in range(P.shape[0]) if i != iq]
    Pq = P[np.ix_(keep, keep)] - np.outer(P[keep, iq], P[iq, keep]) / P[iq, iq]
    Fm = Forms(st, sU); w = Fm.w
    Pi = np.eye(u) - np.outer(w, w) / w.dot(w)
    Xh = Pi @ rng.normal(size=(u, u)) @ Pi
    a = Pi @ rng.normal(size=u); c = rng.normal(size=u); rho = rng.normal(); zeta = rng.normal(size=u)
    Om = np.outer(a, w) + np.outer(w, c)
    def Phiq(X, rho, zeta):
        v = np.concatenate([X.reshape(-1), [rho], zeta]); return v @ Pq @ v
    lhs = Phiq(Xh + Om, rho, zeta)
    FXh = Phiq(Xh, 0, np.zeros(u)); FOm = Phiq(Om, 0, np.zeros(u))
    rhs = FXh + st * Fm.f(a, c) + st * np.sum(Fm.Gcal(a, c) * Xh) + np.sum(Fm.Hfull(rho, zeta) * (Xh + Om)) + Fm.Acal(rho, zeta)
    e1 = abs(FOm - st * Fm.f(a, c)); e2 = abs(lhs - rhs)
    # check Hfull on a random full X
    X = rng.normal(size=(u, u))
    e3 = abs((Phiq(X, rho, zeta) - Phiq(X, 0, 0 * zeta) - Phiq(0 * X, rho, zeta)) / 2 - np.sum(Fm.Hfull(rho, zeta) * X) / 2 * 1)
    e4 = abs(Phiq(0 * X, rho, zeta) - Fm.Acal(rho, zeta))
    mx = max(mx, e1, e2, e3, e4)
    print(u, f'{e1:.2e} {e2:.2e} {e3:.2e} {e4:.2e}')
print('max err', mx)
