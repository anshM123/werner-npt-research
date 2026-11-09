"""Checks of the norm bounds for the map P_a(X) = a Tr(X) I - (X + X^T)/4 (Section 7, "Norm bounds";
results from an independent audit, re-derived and checked here).
 (1) the Hilbert-Schmidt absolute value |P_a| = (a - 1/3) D + (id + T)/4, its Choi matrix (a - 7/12) I + (I + F + Om)/4,
     and the real product identity  I + F + Om = 2 sum_i |ii><ii| + (1/4) sum_{e} |ee><ee|,  e in {+++, +--, -+-, --+};
 (2) max over unit x of ||P_a^{(x)n}(xx^*)||_2^2 = (3a^2 - a + 1/4)^n   (n = 1, 2; optimisation from random starts);
 (3) the four-linear bound |y1^* P_a^{(x)n}(x1 x2^*) y2| <= c_a^n |x1||x2||y1||y2|, c_a = max(a + 1/6, 3/4) (n = 1, 2);
 (4) negative Choi index (9^n - 7^n)/2 and negative-part trace [(7a + 1/2)^n - (9a - 3/2)^n]/2 for 3/4 <= a < 1 (n <= 3);
 (5) no global shift lambda id + mu T makes P_a^{(x)2} completely positive (the two eigenvectors of the audit).
usage: python norm_bounds_check.py"""
import itertools
import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(1)
def P_map(a, X):
    return a * np.trace(X) * np.eye(X.shape[0]) - (X + X.T) / 4
def superop(a):
    """9 x 9 matrix of P_a acting on vec(X) (row-major)"""
    S = np.zeros((9, 9))
    for k in range(9):
        E = np.zeros(9); E[k] = 1
        S[:, k] = P_map(a, E.reshape(3, 3)).reshape(-1)
    return S
def choi(S, d=3):
    C = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            Eij = np.zeros((d, d)); Eij[i, j] = 1
            out = (S @ Eij.reshape(-1)).reshape(d, d)
            C += np.kron(Eij, out)
    return C
Fsw = np.zeros((9, 9))
for i in range(3):
    for j in range(3): Fsw[3 * i + j, 3 * j + i] = 1
om_vec = np.zeros(9); om_vec[[0, 4, 8]] = 1; Om = np.outer(om_vec, om_vec)
# (1)
for a in (0.6, 7 / 12, 0.7, 0.75, 0.9):
    S = superop(a); w, V = np.linalg.eigh((S + S.T) / 2)
    absS = V @ np.diag(np.abs(w)) @ V.T
    Dsup = superop(0.0) * 0; Dm = np.zeros((9, 9))
    for k in range(9):
        E = np.zeros(9); E[k] = 1; X = E.reshape(3, 3); Dm[:, k] = (np.trace(X) * np.eye(3)).reshape(-1)
    Tm = np.zeros((9, 9))
    for k in range(9):
        E = np.zeros(9); E[k] = 1; Tm[:, k] = E.reshape(3, 3).T.reshape(-1)
    pred = (a - 1 / 3) * Dm + (np.eye(9) + Tm) / 4
    Cabs = choi(absS); Cpred = (a - 7 / 12) * np.eye(9) + (np.eye(9) + Fsw + Om) / 4
    print(f"(1) a={a:.4f}: | |P_a| - formula | = {np.abs(absS - pred).max():.1e}, | Choi - formula | = {np.abs(Cabs - Cpred).max():.1e},"
          f" least Choi eigenvalue {np.linalg.eigvalsh(Cabs).min():+.4f} (a - 7/12 = {a - 7/12:+.4f})")
plus = np.ones(3); vecs = []
for sg in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]:
    e = np.array(sg, float); vecs.append(np.kron(e, e))
rhs = 2 * sum(np.outer(np.kron(np.eye(3)[i], np.eye(3)[i]), np.kron(np.eye(3)[i], np.eye(3)[i])) for i in range(3)) \
      + 0.25 * sum(np.outer(v, v) for v in vecs)
print(f"(1) product identity I + F + Om = 2 sum|ii><ii| + (1/4) sum |ee><ee|: error {np.abs(np.eye(9) + Fsw + Om - rhs).max():.1e}")
# helpers for tensor powers
def Pn_apply(a, X, n):
    T = X.reshape((3,) * (2 * n)).astype(complex)
    for s in range(n):
        # apply P_a on copy s
        perm = [s, n + s] + [i for i in range(2 * n) if i not in (s, n + s)]
        Tp = np.transpose(T, perm).reshape(3, 3, -1)
        out = np.empty_like(Tp)
        for k in range(Tp.shape[2]):
            out[:, :, k] = P_map(a, Tp[:, :, k])
        inv = np.argsort(perm)
        T = np.transpose(out.reshape([3, 3] + [3] * (2 * n - 2)), inv)
    return T.reshape(3 ** n, 3 ** n)
# (2) and (3)
for a in (0.6, 0.75, 0.9):
    for n in (1, 2):
        Dn = 3 ** n
        def neg_purity(v):
            x = v[:Dn] + 1j * v[Dn:]; x = x / np.linalg.norm(x)
            Y = Pn_apply(a, np.outer(x, x.conj()), n)
            return -np.vdot(Y, Y).real
        best = max(-minimize(neg_purity, rng.normal(size=2 * Dn), method='BFGS').fun for _ in range(6))
        print(f"(2) a={a}, n={n}: max purity found {best:.8f}, formula (3a^2-a+1/4)^n = {(3*a*a - a + 0.25)**n:.8f}")
        def neg_four(v):
            xs = [v[k * 2 * Dn:(k + 1) * 2 * Dn] for k in range(4)]
            xs = [z[:Dn] + 1j * z[Dn:] for z in xs]; xs = [z / np.linalg.norm(z) for z in xs]
            Y = Pn_apply(a, np.outer(xs[0], xs[1].conj()), n)
            return -abs(xs[2].conj() @ Y @ xs[3])
        bestf = max(-minimize(neg_four, rng.normal(size=8 * Dn), method='BFGS').fun for _ in range(6))
        ca = max(a + 1 / 6, 0.75)
        print(f"(3) a={a}, n={n}: max four-linear value found {bestf:.6f}, bound c_a^n = {ca ** n:.6f}")
# (4)
a = 0.8
C1 = choi(superop(a)); ev1 = np.linalg.eigvalsh(C1)
print(f"(4) one-copy Choi eigenvalues at a={a}: {np.round(np.sort(ev1), 4).tolist()} (expect a-1 once, a-1/4 x5, a+1/4 x3)")
for n in (1, 2, 3):
    ev = np.array([1.0])
    for _ in range(n): ev = np.kron(ev, ev1)
    neg = ev[ev < -1e-12]
    print(f"(4) n={n}: negative Choi eigenvalues {len(neg)} (formula {(9**n - 7**n)//2}), their total {-neg.sum():.6f}"
          f" (formula {((7*a + 0.5)**n - (9*a - 1.5)**n)/2:.6f})")
# (5) global shifts at n = 2.  Global Choi matrix of P_a (x) P_a, index (in1 in2, out1 out2), built directly from the map
def choi_global(a):
    C = np.zeros((81, 81), dtype=complex)
    for I in range(9):
        for J in range(9):
            E = np.zeros((9, 9)); E[I, J] = 1
            C += np.kron(E, Pn_apply(a, E, 2))
    return C.real
C2 = choi_global(a)
Fbig = np.zeros((81, 81))
for i in range(9):
    for j in range(9): Fbig[9 * i + j, 9 * j + i] = 1      # Choi matrix of the global transpose
omb = np.zeros(81); omb[[10 * i for i in range(9)]] = 1    # Choi vector of the global identity
o = om_vec / np.sqrt(3)
symt = (np.kron(np.eye(3)[0], np.eye(3)[0]) - np.kron(np.eye(3)[1], np.eye(3)[1])) / np.sqrt(2)      # symmetric traceless
anti = (np.kron(np.eye(3)[0], np.eye(3)[1]) - np.kron(np.eye(3)[1], np.eye(3)[0])) / np.sqrt(2)
def to_global(v1, v2):
    """per-copy Choi vectors v1 (in1 out1), v2 (in2 out2) -> global index (in1 in2, out1 out2)"""
    T = np.kron(v1, v2).reshape(3, 3, 3, 3)          # (in1, out1, in2, out2)
    return np.transpose(T, (0, 2, 1, 3)).reshape(-1)
C1 = choi(superop(a))
Pperm = np.zeros((81, 81))
for k in range(81):
    e = np.zeros(81); e[k] = 1; Pperm[:, k] = np.transpose(e.reshape(3, 3, 3, 3), (0, 2, 1, 3)).reshape(-1)
print(f"(5) global Choi of P_a (x) P_a = reordered kron(Choi, Choi): error {np.abs(C2 - Pperm @ np.kron(C1, C1) @ Pperm.T).max():.1e}")
vals = {}
for name, v in [("omega x symmetric", to_global(o, symt)), ("omega x antisymmetric", to_global(o, anti))]:
    vals[name] = (v @ C2 @ v, v @ Fbig @ v, v @ omb)
    print(f"(5) {name}: Choi value {vals[name][0]:+.5f}, <v|F|v> = {vals[name][1]:+.2f}, <v|global Omega> = {vals[name][2]:+.1e}")
print(f"(5) expected (a-1)(a-1/4) = {(a-1)*(a-0.25):+.5f} and (a-1)(a+1/4) = {(a-1)*(a+0.25):+.5f}; adding lambda id + mu T"
      f" needs mu >= {-vals['omega x symmetric'][0]:.4f} and mu <= {vals['omega x antisymmetric'][0]:.4f}: impossible")
