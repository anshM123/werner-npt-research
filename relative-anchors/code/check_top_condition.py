"""Checks for Theorem K (README, Section 3c): the top group condition (R_inf)_k for every k and all weights.

Torus blocks of (R_inf)_k: for each set T of copies whose column letter is 2, the type I block on
P_r = (1@r, y_T), r in [k], and Z_t = (2@t, y_T with t -> 1), t in T; and the one-state type II blocks (2@t, y_B), t in B.
(a) The closed form (P, C, B; type II scalars) against a direct construction of Y from M, V, D (k = 2..5, every T).
(b) Exact identities at rational weights: C = P L with L_t = -a e_t/(2 a_t), and F = B - L^T P L has
    F_tt = 1 - a(a - s_t)/(4 a_t^2), F_tu = -w_t w_u s_t s_u/(16 a_t a_u a_tu); diagonal dominance margin of F >= 11/16.
(c) The gap: on the complement of (w, 0) every type I block is >= 11/48, and every type II scalar is >= (1 - max s/|s|)/2.
(d) Proposition M: the 1 x 1 condition (R_inf)_{-1} equals 2^{-k}/(1 + gamma/2), gamma = sum_r s_r/(2(|s| - s_r)).
Uses check_bath_reduction.group_ops from single-excitation-range/code for (a)."""
import os, sys, itertools
import numpy as np
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', 'single-excitation-range', 'code'))
from check_bath_reduction import group_ops


def blocks(w, T, F=float):
    k = len(w); s = [x * x for x in w]; S = sum(s); ST = sum(s[t] for t in T); T = list(T)
    N = [F(1) if r in T else F(1) / 2 for r in range(k)]; v = [w[r] * N[r] for r in range(k)]
    a = (S + ST) / 2; at = {t: a - s[t] / 2 for t in T}
    P = [[(N[r] if r == q else F(0)) - v[r] * v[q] / a for q in range(k)] for r in range(k)]
    C = [[(a / (2 * at[t])) * (-(F(1) if r == t else F(0)) + w[t] * v[r] / a) for t in T] for r in range(k)]
    B = [[F(1) if t == u else -w[t] * w[u] / (4 * (a - (s[t] + s[u]) / 2)) for u in T] for t in T]
    return P, C, B, a, at, s


def assemble(P, C, B, k, m):
    return np.block([[np.array(P, float), np.array(C, float).reshape(k, m)],
                     [np.array(C, float).reshape(k, m).T, np.array(B, float).reshape(m, m)]])


def direct(w, T):
    k = len(w); pG = {}
    for t in range(k):
        x = [0] * k; x[t] = 1; pG[tuple(x)] = w[t]
    _, Aent, Bent = group_ops(pG, k)
    yT = lambda V: tuple(2 if t in V else 1 for t in range(k))
    st = [(tuple(1 if t == r else 0 for t in range(k)), yT(T)) for r in range(k)]
    st += [(tuple(2 if x == t else 0 for x in range(k)), yT(set(T) - {t})) for t in T]
    full = np.array([[Aent(a, yp, b, y) - Bent(a, yp, b, y) for (b, y) in st] for (a, yp) in st])
    two = {t: Aent(tuple(2 if x == t else 0 for x in range(k)), yT(T), tuple(2 if x == t else 0 for x in range(k)), yT(T))
           - Bent(tuple(2 if x == t else 0 for x in range(k)), yT(T), tuple(2 if x == t else 0 for x in range(k)), yT(T))
           for t in T}
    return full, two


ok = True
rng = np.random.default_rng(11)
for k in (2, 3, 4, 5):
    w = list(rng.random(k) + 0.1); err = 0; nst = 0
    for m in range(k + 1):
        for T in itertools.combinations(range(k), m):
            P, C, B, a, at, s = blocks(w, T)
            full, two = direct(w, T)
            err = max(err, np.abs(full - assemble(P, C, B, k, m)).max())
            for t in T:
                err = max(err, abs(two[t] - (a - s[t]) / (2 * a - s[t])))
            nst += k + 2 * m
    print(f"(a) k={k}: closed form vs direct construction, max error {err:.1e}; states covered {nst} = 2k 2^k = {2 * k * 2 ** k}")
    ok &= err < 1e-12 and nst == 2 * k * 2 ** k
for trial in range(40):
    k = int(rng.integers(2, 7)); w = [Fr(int(rng.integers(1, 30)), int(rng.integers(1, 30))) for _ in range(k)]
    m = int(rng.integers(1, k + 1)); T = sorted(rng.choice(k, m, replace=False).tolist())
    P, C, B, a, at, s = blocks(w, T, Fr)
    L = [[(-a / (2 * at[t]) if r == t else Fr(0)) for t in T] for r in range(k)]
    PL = [[sum(P[r][q] * L[q][j] for q in range(k)) for j in range(m)] for r in range(k)]
    ok &= all(PL[r][j] == C[r][j] for r in range(k) for j in range(m))
    LPL = [[sum(L[q][i] * PL[q][j] for q in range(k)) for j in range(m)] for i in range(m)]
    for i, t in enumerate(T):
        row = [B[i][j] - LPL[i][j] for j in range(m)]
        for j, u in enumerate(T):
            if t == u: ok &= row[j] == 1 - a * (a - s[t]) / (4 * at[t] ** 2)
            else: ok &= row[j] == -w[t] * w[u] * s[t] * s[u] / (16 * at[t] * at[u] * (a - (s[t] + s[u]) / 2))
        ok &= row[i] - sum(abs(row[j]) for j in range(m) if j != i) >= Fr(11, 16)
print(f"(b) exact identities (C = PL, entries of F, dominance margin >= 11/16) at 40 rational weight vectors: {'PASS' if ok else 'FAIL'}")
worst = 9.0; worst2 = 9.0
for trial in range(4000):
    k = int(rng.integers(2, 10)); w = rng.random(k) ** rng.choice([1, 4, 10]) + 1e-4
    if trial % 9 == 0: w = np.array([1.0] + [1e-3] * (k - 1))
    s = w ** 2; nu2 = 0.5 * (1 - s.max() / s.sum())
    for m in range(k + 1):
        if trial > 30 and rng.random() > 0.3: continue
        T = sorted(rng.choice(k, m, replace=False).tolist())
        P, C, B, a, at, _ = blocks(list(w), T)
        Y = assemble(P, C, B, k, m)
        e = np.zeros(k + m); e[:k] = w / np.linalg.norm(w); Pi = np.eye(k + m) - np.outer(e, e)
        ok &= abs(e @ Y @ e) < 1e-12
        worst = min(worst, np.linalg.eigvalsh(Pi @ Y @ Pi + 10 * np.outer(e, e))[0])
        for t in T: worst2 = min(worst2, (a - s[t]) / (2 * a - s[t]) - nu2)
print(f"(c) smallest type I eigenvalue off (w, 0): {worst:.4f} (>= 11/48 = {11/48:.4f}); type II minus (1 - max s/|s|)/2 >= {worst2:.2e}")
ok &= worst >= 11 / 48 - 1e-12 and worst2 >= -1e-15
err = 0
for k in (2, 3, 4, 5, 6):
    for trial in range(3):
        w = rng.random(k) ** (1 + 3 * trial) + 0.01; s = w ** 2; pG = {}
        for t in range(k):
            x = [0] * k; x[t] = 1; pG[tuple(x)] = w[t]
        _, Aent, Bent = group_ops(pG, k); vac = (0,) * k
        gam = sum(s[r] / (2 * (s.sum() - s[r])) for r in range(k))
        err = max(err, abs(Aent(vac, vac, vac, vac) - Bent(vac, vac, vac, vac) - 2.0 ** (-k) / (1 + gam / 2)))
print(f"(d) Proposition M: (R_inf)_(-1) = 2^(-k)/(1 + gamma/2), k = 2..6, max error {err:.1e}")
ok &= err < 1e-13
print("PASS" if ok else "FAIL")
