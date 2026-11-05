"""Exploration for the (R_inf)_2 blocks with |T| = 3 (no letter 1): closed form against Lemma B (j1lib.entry works for
every j) and ylib, and the vertex relaxation (one Stieltjes sum exact, joint concavity) in floating point."""
import itertools
import numpy as np
from fractions import Fraction as Fr
import flint
from j1lib import J1, onehot, mu, to_ylib_state
from ylib import Builder, point_ring


def states_j2T3(k, T):
    t1, t2, t3 = T
    st = []
    for i, t in enumerate(T):
        rest = [u for u in T if u != t]
        st.append(((2, t), onehot(k, {rest[0]: 2, rest[1]: 2})))
    st.append((None, onehot(k, {t1: 2, t2: 2, t3: 2})))
    return st


def closed(x, T):
    k = len(x)
    U = [u for u in range(k) if u not in T]
    xT = sum(x[t] for t in T)
    M = np.zeros((4, 4))
    for i in range(3):
        M[i, i] = 2
    for i, j in itertools.combinations(range(3), 2):
        kk = 3 - i - j
        tk = T[kk]
        g = lambda u: x[u] / (2 - 2 * x[u] - x[tk])
        ai = 1 + g(T[i]); aj = 1 + g(T[j])
        Gam = 1 + 0.5 * sum(g(u) for u in range(k) if u != tk)
        M[i, j] = M[j, i] = ai * aj / (2 * Gam)
    for i in range(3):
        rest = [T[l] for l in range(3) if l != i]
        sh = x[rest[0]] + x[rest[1]]
        g = lambda u: x[u] / (2 - 2 * x[u] - sh)
        api = 1 + g(T[i])
        Gp = 1 + 0.5 * (g(T[i]) + sum(g(u) for u in U))
        M[i, 3] = M[3, i] = -api / Gp
    GT = 1 + 0.5 * sum(x[u] / (2 - 2 * x[u] - xT) for u in U)
    M[3, 3] = 2 / GT
    return M


ok = True
rng = np.random.default_rng(1)
for k in (4, 5, 6):
    w = [Fr(int(a), int(b)) for a, b in zip(rng.integers(1, 9, k), rng.integers(1, 9, k))]
    J = J1(w)
    T = (0, 1, 2) if k > 4 else (1, 2, 3)
    st = states_j2T3(k, T)
    Mx = [[J.entry(a, b) for b in st] for a in st]
    S = sum(J.s)
    x = [float(s / S) for s in J.s]
    C = closed(x, T)
    err = np.abs(np.array(Mx, dtype=float) - C).max()
    Bd = Builder(point_ring([flint.fmpq(v.numerator ** 2, v.denominator ** 2) for v in w]))
    ys = [to_ylib_state(z, k) for z in st]
    same = all(Fr(int(Bd.Y(ys[i], ys[j]).p), int(Bd.Y(ys[i], ys[j]).q)) == Mx[i][j] / 2 ** (k - 2) for i in range(4) for j in range(4))
    print(f'k={k}: closed form error {err:.2e}; Lemma B entries (units 2^(k-2)) = ylib: {same}')
    ok &= err < 1e-12 and same


# vertex relaxation
def h(d, z):
    return z / (1 + d * z)


def vertex_min(xt, e):
    """xt = (x_t1, x_t2, x_t3), e = mass of U (sum 1).  Returns min over lam grid and the 64 vertices of
    f = 2 - G_T a^T C^{-1} a, and whether C > 0 at the 8 c-vertices."""
    xT = sum(xt)
    Z0, Z1 = 1 / (2 - xT), 1 / xT
    out = 1e9; cpd = True
    for lam in np.linspace(0, 1, 21):
        Z = Z0 + lam * (Z1 - Z0)
        GT = 1 + e * Z / 2
        cr = []
        for i, j in itertools.combinations(range(3), 2):
            kk = 3 - i - j
            d = xT - xt[kk]
            Flo = e * ((1 - lam) * h(d, Z0) + lam * h(d, Z1)); Fhi = e * h(d, Z)
            g = lambda u: u / (2 - 2 * u - xt[kk])
            ai = 1 + g(xt[i]); aj = 1 + g(xt[j])
            G0 = 1 + 0.5 * (g(xt[i]) + g(xt[j]))
            cr.append(sorted([ai * aj / (2 * (G0 + F / 2)) for F in (Flo, Fhi)]))
        ar = []
        for i in range(3):
            d = xt[i]
            sh = xT - xt[i]
            Flo = e * ((1 - lam) * h(d, Z0) + lam * h(d, Z1)); Fhi = e * h(d, Z)
            gi = xt[i] / (2 - 2 * xt[i] - sh)
            ar.append(sorted([(1 + gi) / (1 + 0.5 * (gi + F)) for F in (Flo, Fhi)]))
        for cv in itertools.product(*cr):
            Cm = np.array([[2, cv[0], cv[1]], [cv[0], 2, cv[2]], [cv[1], cv[2], 2]])
            if np.linalg.eigvalsh(Cm)[0] <= 0:
                cpd = False
                continue
            Ci = np.linalg.inv(Cm)
            for av in itertools.product(*ar):
                a = np.array(av)
                out = min(out, 2 - GT * a @ Ci @ a)
    return out, cpd


worst = 1e9; allpd = True
for trial in range(1500):
    p = rng.dirichlet([1, 1, 1, 1] if trial % 2 else [0.3, 0.3, 0.3, 0.3])
    if p[3] < 1e-9:
        continue
    v, cpd = vertex_min(p[:3], p[3])
    allpd &= cpd
    r = v / (1 - p[:3].max())
    if r < worst:
        worst = r; info = p
print('vertex relaxation: min f/(1 - max x_t) =', worst, info, '; C > 0 at all c-vertices:', allpd)
# true configurations for comparison
worst2 = 1e9
for trial in range(1500):
    k = int(rng.integers(4, 9))
    x = rng.dirichlet(np.ones(k) * rng.choice([0.3, 1, 3]))
    M = closed(list(x), (0, 1, 2))
    Cm = M[:3, :3]; a = -M[:3, 3]; GT = 2 / M[3, 3]
    f = 2 - GT * a @ np.linalg.solve(Cm, a)
    worst2 = min(worst2, f / (1 - x[:3].max()))
print('true blocks: min f/(1 - max x_t) =', worst2)
print('PASS' if ok else 'FAIL')
