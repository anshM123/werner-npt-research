"""Checks for the exact reduction of the type I blocks of (R_inf)_{k-1} (README, end of Section 3) and for the constant
5/24 in step 5 of Theorem T.
(a) The principal submatrix on P_{t,h} (row letter 1 at t in T, column: 2 on T, 1 on U with a hole at h in U) equals
    K = 1/2 I - 1/2 w_T w_T^T (x) E^{-1}, E = S_T I + E_0, E_0 = diag(s_h + S_U/2) + 1/2 w_U w_U^T  (direct construction,
    k = 3, 4, 5, every |T|); K^{-1} = 2I + 2 w_T w_T^T (x) E_0^{-1}; K >= S_U/(2(2 S_T + S_U)); the block has
    u(k + m) + m + 1 states.
(b) Theorem T, step 5: the form F has lambda_min >= 5/24 (random and extreme weights), and the scalar lemma
    (1 + 2x)(1 + 2y)(2 + rho_a + rho_b - 12 rho_a rho_b) = 2 + 6t - 32p >= 2(1 - t)(1 + 4t) for p <= t^2/4 (sympy)."""
import os, sys, itertools
import numpy as np, sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', 'single-excitation-range', 'code'))
from check_bath_reduction import group_ops

ok = True
rng = np.random.default_rng(2)
for k in (3, 4, 5):
    w = rng.random(k) + 0.2; s = w ** 2; pG = {}
    for t in range(k):
        x = [0] * k; x[t] = 1; pG[tuple(x)] = w[t]
    _, Aent, Bent = group_ops(pG, k)
    for m in range(1, k):
        T = list(range(m)); U = list(range(m, k)); u = len(U)
        col = lambda h: tuple(2 if c in T else (0 if c == h else 1) for c in range(k))
        st = [(tuple(1 if c == t else 0 for c in range(k)), col(h)) for t in T for h in U]
        Y = np.array([[Aent(a, yp, b, y) - Bent(a, yp, b, y) for (b, y) in st] for (a, yp) in st])
        ST, SU, wT, wU = s[T].sum(), s[U].sum(), w[T], w[U]
        E0 = np.diag(s[U] + SU / 2) + 0.5 * np.outer(wU, wU); E = ST * np.eye(u) + E0
        K = 0.5 * np.eye(m * u) - 0.5 * np.kron(np.outer(wT, wT), np.linalg.inv(E))
        Kinv = 2 * np.eye(m * u) + 2 * np.kron(np.outer(wT, wT), np.linalg.inv(E0))
        # size of the whole torus block: P (k u), R (m), Z (m u), Q (1)
        yT = tuple(2 if c in T else 1 for c in range(k))
        size = 0
        for a in [tuple(1 if c == r else 0 for c in range(k)) for r in range(k)] + [tuple(2 if c == r else 0 for c in range(k)) for r in range(k)] + [(0,) * k]:
            for y in itertools.product(range(3), repeat=k):
                nz = sum(v != 0 for v in y)
                if nz != (k - 1 if any(a) else k): continue
                twos = [int(a[c] == 2) + int(y[c] == 2) for c in range(k)]
                ones = sum(v == 1 for v in a) + sum(v == 1 for v in y)
                if twos == [1 if c in T else 0 for c in range(k)] and ones == u:
                    size += 1
        e1 = np.abs(Y - K).max(); e2 = np.abs(K @ Kinv - np.eye(m * u)).max(); lm = np.linalg.eigvalsh(K)[0]
        lb = SU / (2 * (2 * ST + SU))
        good = e1 < 1e-12 and e2 < 1e-12 and lm >= lb - 1e-12 and size == u * (k + m) + m + 1
        ok &= good
        print(f"(a) k={k} |T|={m}: |Y - K| = {e1:.1e}, |K K^-1 - I| = {e2:.1e}, lambda_min(K) = {lm:.4f} >= {lb:.4f}, block size {size} = u(k+m)+m+1: {good}")
worst = 9.0
for trial in range(3000):
    k = int(rng.integers(2, 8)); s = rng.random(k) ** rng.choice([1, 3, 8]) + 1e-6
    if trial % 10 == 0: s = np.array([1, 1] + [1e-6] * (k - 2))
    S = s.sum(); rho = s / (s + S / 2); G = rho.sum(); Q = np.zeros((k * k, k * k))
    for a in range(k):
        for b in range(k):
            Q[a * k + b, a * k + b] += 0.25 if a != b else 0.75 - rho[a] / 2
            Q[a * k + b, b * k + a] += rho[a] * rho[b] / (2 * (2 + G))
    worst = min(worst, np.linalg.eigvalsh((Q + Q.T) / 2)[0])
x, y, t = sp.symbols('x y t', positive=True)
ra, rb = 2 * x / (1 + 2 * x), 2 * y / (1 + 2 * y)
lhs = sp.expand(sp.cancel((2 + ra + rb - 12 * ra * rb) * (1 + 2 * x) * (1 + 2 * y)))
lem = sp.expand(lhs - (2 + 6 * (x + y) - 32 * x * y)) == 0 and sp.expand(2 + 6 * t - 8 * t ** 2 - 2 * (1 - t) * (1 + 4 * t)) == 0
print(f"(b) Theorem T form: smallest lambda_min = {worst:.6f} (claim >= 5/24 = {5/24:.6f}); scalar lemma identities: {lem}")
ok &= worst >= 5 / 24 - 1e-9 and lem
print("PASS" if ok else "FAIL")
