"""Checks for Theorem B0 (README, Section 3b): the bottom group condition (R_inf)_0 holds for every k and all weights.

(a) Closed form of (R_inf)_0 (scaled by 2^{k-1}) against a direct construction (k = 2..5).
(b) Letter-2 blocks: the 2 x 2 form [[1, -(1-kappa_r)/2], [-(1-kappa_r)/2, eta_r]] with
    eta_r = 1/(1 + gamma_r/2), 1 - kappa_r = (1 + g_r)/(1 + gamma/2), against the closed form, and the inequality
    (1 + gamma/2)^2 >= (1/4)(1 + g_r)^2 (1 + gamma_r/2) at random and extreme weights.
(c) Letter-1 block: PSD <=> (C0); the block's sign agrees with the (C0) ratio at random weights.
(d) The bound chains of Regimes A and B hold at random (x, zeta) with sum_r x_r zeta_r = 0.
(e) The two reduced two-parameter inequalities, in interval arithmetic (c0_interval.py).
Uses check_bath_reduction.group_ops from single-excitation-range/code for (a)."""
import os, sys, itertools
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', 'single-excitation-range', 'code'))
from check_bath_reduction import group_ops

def strings(k, c): return [s for s in itertools.product(range(3), repeat=k) if sum(v != 0 for v in s) == c]

def closed_form(w):
    k = len(w); s = w ** 2; S = s.sum()
    H1 = strings(k, 1); H2 = strings(k, 2); i1 = {x: i for i, x in enumerate(H1)}; i2 = {x: i for i, x in enumerate(H2)}
    Jp = np.zeros((len(H2), len(H1)))
    for y in H1:
        for t in range(k):
            if y[t] == 0:
                z = list(y); z[t] = 1; Jp[i2[tuple(z)], i1[y]] += w[t]
    def Cd(H): return np.array([sum(2 * s[t] if x[t] == 0 else (s[t] if x[t] == 2 else 0) for t in range(k)) for x in H])
    wv = np.array([w[y.index(1)] if 1 in y else 0.0 for y in H1])
    E1 = np.diag(Cd(H1)) + 0.5 * np.outer(wv, wv); E2 = np.diag(Cd(H2)) + 0.5 * Jp @ Jp.T
    u = np.linalg.solve(E1, wv)
    Y11 = np.eye(len(H1)) - np.outer(wv, wv) / S
    Y00 = np.eye(len(H1)) - 0.5 * Jp.T @ np.linalg.solve(E2, Jp)
    K = np.zeros((len(H1), len(H1)))
    for ia, a in enumerate(H1):
        r = [t for t in range(k) if a[t]][0]; l = a[r]
        for iy, y in enumerate(H1):
            rp = [t for t in range(k) if y[t]][0]; lp = y[rp]
            val = 0.5 * (a == y) * (u @ wv)
            if l == 1: val += w[r] * (2 * (rp != r) + (rp == r and lp == 2)) * u[iy]
            if l == 2 and rp == r and lp == 2: val -= w[r] * u[i1[tuple(1 if t == r else 0 for t in range(k))]]
            K[ia, iy] = val
    Y10 = -0.5 * (np.eye(len(H1)) - K)
    return H1, np.block([[Y11, Y10], [Y10.T, Y00]])

def direct(w):
    k = len(w); pG = {}
    for t in range(k):
        x = [0] * k; x[t] = 1; pG[tuple(x)] = w[t]
    _, Aent, Bent = group_ops(pG, k)
    H1 = strings(k, 1); vac = (0,) * k; st = [(a, vac) for a in H1] + [(vac, y) for y in H1]
    return 2 ** (k - 1) * np.array([[Aent(a, yp, b, y) - Bent(a, yp, b, y) for (b, y) in st] for (a, yp) in st])

def gam(x): return np.sum(x / (2 * (1 - x)))

def c0_parts(x):
    k = len(x); lam = 1 / (1 - x); Lm = np.diag(x * lam ** 2)
    for r, t in itertools.combinations(range(k), 2):
        v = np.zeros(k); v[r] = lam[r]; v[t] = lam[t]
        Lm += 0.25 * x[r] * x[t] / (1 - x[r] - x[t]) * np.outer(v, v)
    return lam, Lm

rng = np.random.default_rng(2026); ok = True
# (a)
for k in (2, 3, 4, 5):
    w = rng.random(k) + 0.1; _, Yc = closed_form(w); err = np.abs(Yc - direct(w)).max()
    print(f"(a) k={k}: closed form vs direct, max error {err:.1e}"); ok &= err < 1e-10
# (b)
worst = 0
for trial in range(3000):
    k = rng.integers(2, 9); x = rng.dirichlet(np.ones(k) * rng.choice([0.2, 1, 5]))
    if trial % 10 == 0: x = np.array([1 - 1e-6] + [1e-6 / (k - 1)] * (k - 1))
    x = np.maximum(x, 1e-12); x /= x.sum(); g = gam(x)
    for r in range(k):
        gr = x[r] / (2 * (1 - x[r])); grr = sum(x[t] / (2 - x[r] - 2 * x[t]) for t in range(k) if t != r)
        worst = max(worst, 0.25 * (1 + gr) ** 2 * (1 + grr / 2) / (1 + g / 2) ** 2)
    if trial < 20 and 2 <= k <= 4:
        H1, Yc = closed_form(np.sqrt(x)); m = len(H1)
        for r in range(k):
            i_r = H1.index(tuple(2 if t == r else 0 for t in range(k)))
            gr = x[r] / (2 * (1 - x[r])); grr = sum(x[t] / (2 - x[r] - 2 * x[t]) for t in range(k) if t != r)
            blk = Yc[np.ix_([i_r, m + i_r], [i_r, m + i_r])]
            ok &= abs(blk[1, 1] - 1 / (1 + grr / 2)) < 1e-10 and abs(-2 * blk[0, 1] - (1 + gr) / (1 + g / 2)) < 1e-10 and abs(blk[0, 0] - 1) < 1e-12
print(f"(b) letter-2 blocks: formula matches the closed form; largest ratio (1/4)(1+g_r)^2(1+gamma_r/2)/(1+gamma/2)^2 = {worst:.6f} (<= 1)")
ok &= worst <= 1 + 1e-9
# (c) and (d)
viol_c = 0; viol_d = 0
for trial in range(2000):
    k = rng.integers(3, 8); x = rng.dirichlet(np.ones(k) * rng.choice([0.2, 1, 5]))
    if trial % 7 == 0: x = np.array([0.5, 0.5 - 1e-4] + [1e-4 / (k - 2)] * (k - 2))
    x = np.maximum(x, 1e-12); x /= x.sum(); xs = np.sort(x)[::-1]; g = gam(x)
    lam, Lm = c0_parts(x)
    # (c): max of L/S on the constraint hyperplane vs PSD of the letter-1 block
    Q, _ = np.linalg.qr(np.column_stack([x] + [rng.normal(size=k) for _ in range(k - 1)])); B = Q[:, 1:]
    Sr = B.T @ np.diag(x) @ B; ws, Vs = np.linalg.eigh(Sr); Sh = Vs @ np.diag(ws ** -0.5) @ Vs.T
    ratio = np.linalg.eigvalsh(Sh @ B.T @ Lm @ B @ Sh)[-1] / (2 + g) ** 2
    if trial < 200 and k <= 5:
        H1, Yc = closed_form(np.sqrt(x)); m = len(H1)
        l1 = [i for i, a in enumerate(H1) if 2 not in a]; idx = l1 + [m + i for i in l1]
        psd = np.linalg.eigvalsh(Yc[np.ix_(idx, idx)])[0] > -1e-10
        viol_c += (psd != (ratio <= 1 + 1e-9))
    # (d): bound chains at random zeta on the hyperplane
    for _ in range(5):
        z = B @ rng.normal(size=k - 1); S = np.sum(x * z ** 2); L = z @ Lm @ z
        q = x * z ** 2 / S; order = np.argsort(-x); a, b = xs[0], xs[1]
        if a + b <= 0.75:
            bound = np.sum(q * lam * (lam + 2))
        else:
            l1_, l2_ = lam[order[0]], lam[order[1]]; eps = 1 - a - b
            c = np.empty(k)
            for r in range(k):
                if r == order[0]: c[r] = lam[r] ** 2 + l2_ ** 2
                elif r == order[1]: c[r] = lam[r] ** 2 + lam[r] * min(l1_, 1 / (4 * x[r])) + 4 / 3 * lam[r] ** 2 * eps
                else: c[r] = lam[r] ** 2 + lam[r] * min(l1_, 1 / (4 * x[r])) + 4 / 3 * lam[r] ** 2 * (1 - a - x[r])
            bound = np.sum(c * q) + 0.5 * l1_ * l2_
        viol_d += L / S > bound + 1e-9
print(f"(c) PSD of the letter-1 block <=> (C0): disagreements {viol_c}"); ok &= viol_c == 0
print(f"(d) bound chains of Regimes A and B: violations {viol_d}"); ok &= viol_d == 0
# (e)
import c0_interval
from mpmath import mpf
okA, infoA = c0_interval.run(c0_interval.checkA, ((mpf(0), mpf(3) / 4), (mpf(0), mpf(3) / 8)), c0_interval.insideA, maxdepth=14)
okB, infoB = c0_interval.run(c0_interval.checkB, ((mpf(3) / 4, mpf(1)), (mpf(0), mpf(1) / 2)), c0_interval.insideB, maxdepth=14)
print(f"(e) interval arithmetic: Regime A {'PASS' if okA else 'FAIL'} {infoA}; Regime B (x1 <= 0.99) {'PASS' if okB else 'FAIL'} {infoB}")
ok &= okA and okB
print("PASS" if ok else "FAIL")
