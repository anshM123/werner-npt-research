"""Counterexamples to Conjecture SC and SC-infinity in the tangent class.

Range u_z = sum_t |z at copy t> (z = 1, 2): U_t = I at every copy and a = 0 (the W-type single-excitation plane).
Co-range: integer (n = 3) or Gaussian-integer (n = 4, 5, 6) vectors g_1, g_2, and C = u_1 g_1^T + u_2 g_2^T.
At copy 0 the pencil Q_0 - lambda N_0 is positive definite at lambda = -1/2, so it is PSD for some lambda in [-1/2, 0)
(against Conjecture SC, n >= 4) and there is no H with N_0(H) = 0, Q_0(H) < 0 (against SC-infinity).
n = 3 is checked in exact rational arithmetic (sympy, Sylvester's criterion on an unnormalised basis of Herm_3);
n = 4, 5, 6 in floating point, with the forms checked against direct evaluation of q(GCG) and ||GCG||^2.
Every example has rank two and is thick at every copy (range and co-range).  Needs sympy."""
import numpy as np, itertools, sympy as sp
from scipy.linalg import sqrtm
from common import squeeze_forms, q, q_direct, tangent_range, Bt

rng = np.random.default_rng(0)


def local_ratios(C, n):
    D = 3 ** n; out = []
    for s in range(n):
        for Y in (C @ C.conj().T, C.conj().T @ C):
            T = np.moveaxis(Y.reshape([3] * (2 * n)), [s, n + s], [0, n]).reshape(3, D // 3, 3, D // 3)
            ev = np.linalg.eigvalsh(np.einsum('aibi->ab', T)); out.append(ev[0] / ev[-1])
    return min(out)


def embed(G, s, n):
    ops = [np.eye(3)] * n; ops[s] = G; out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def float_report(C, n):
    C = C / np.linalg.norm(C); K, Kn = squeeze_forms(C, n, 0); err = 0
    for _ in range(3):
        A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); M = A @ A.conj().T; G = embed(sqrtm(M), 0, n)
        Y = G @ C @ G; c = np.array([np.trace(B.conj().T @ M).real for B in Bt])
        err = max(err, abs(c @ K @ c - q_direct(Y, n)) / abs(q_direct(Y, n)), abs(c @ Kn @ c - np.linalg.norm(Y) ** 2) / np.linalg.norm(Y) ** 2)
    lams = np.linspace(-4, 1, 5001); vals = np.array([np.linalg.eigvalsh(K - l * Kn)[0] for l in lams]); pd = lams[vals > 0]
    margin = np.linalg.eigvalsh(K + Kn / 2)[0] / np.linalg.norm(Kn, 2)
    others = [(lambda K_, Kn_: max(np.linalg.eigvalsh(K_ - l * Kn_)[0] for l in lams[::10]) / np.linalg.norm(Kn_, 2))(*squeeze_forms(C, n, s)) for s in range(1, n)] if n <= 4 else []
    print(f"n={n}: rank {np.linalg.matrix_rank(C, 1e-10)}, smallest local-support ratio {local_ratios(C, n):.3f}, R = {q(C, n):.4f}; "
          f"forms vs direct {err:.1e}; copy 0: lambda_min(Q + N/2)/||N|| = {margin:+.3e}, PD for lambda in [{pd.min():.2f}, {pd.max():.2f}]"
          + (f"; other copies: largest margin {max(others):+.3f}" if others else ""), flush=True)
    return margin > 1e-5 and err < 1e-10 and np.linalg.matrix_rank(C, 1e-10) == 2 and local_ratios(C, n) > 0.05


ok = True
# ---- n = 3: explicit integer example, exact ----
n = 3; D = 27
def ket(s):
    v = np.zeros(D, dtype=np.int64); v[9 * s[0] + 3 * s[1] + s[2]] = 1; return v
u = tangent_range([np.eye(2)] * n, np.zeros(2), n).real.round().astype(np.int64)
g1 = -ket((0,0,0)) - 3*ket((0,0,1)) - 3*ket((0,1,0)) + ket((0,2,0)) + 3*ket((1,0,0)) - ket((1,0,1)) - ket((1,1,2)) - ket((1,2,1)) + 2*ket((2,1,1))
g2 = ket((0,0,0)) - 3*ket((0,0,2)) - 3*ket((0,2,0)) - 2*ket((1,2,2)) + 3*ket((2,0,0)) + ket((2,0,2)) + ket((2,1,2)) + ket((2,2,1))
C3 = np.outer(u[0], g1) + np.outer(u[1], g2)
ok &= float_report(C3.astype(float), 3)
Cs = sp.Matrix(C3.tolist())


def ptrace(Mx, J, k):
    rest = [t for t in range(k) if t not in J]; dr = 3 ** len(rest); out = sp.zeros(dr, dr)
    for i in range(3 ** k):
        si = np.unravel_index(i, [3] * k)
        for j in range(3 ** k):
            if Mx[i, j] == 0:
                continue
            sj = np.unravel_index(j, [3] * k)
            if any(si[t] != sj[t] for t in J):
                continue
            ii = int(np.ravel_multi_index([si[t] for t in rest], [3] * len(rest))) if rest else 0
            jj = int(np.ravel_multi_index([sj[t] for t in rest], [3] * len(rest))) if rest else 0
            out[ii, jj] += Mx[i, j]
    return out, rest


def at_copy(Hm, s, k):
    mats = [sp.eye(3)] * k; mats[s] = Hm; out = mats[0]
    for m_ in mats[1:]:
        out = sp.kronecker_product(out, m_)
    return out


basis = []
for i in range(3):
    E = sp.zeros(3, 3); E[i, i] = 1; basis.append(E)
for i in range(3):
    for j in range(i + 1, 3):
        E = sp.zeros(3, 3); E[i, j] = E[j, i] = 1; basis.append(E)
        E = sp.zeros(3, 3); E[i, j] = -sp.I; E[j, i] = sp.I; basis.append(E)
pieces = []
for k in range(3):
    for Jp in itertools.combinations([1, 2], k):
        Y, rest = ptrace(Cs, list(Jp), 3); pieces.append((sp.Rational(-1, 2) ** k, Y, rest.index(0), len(rest), list(Jp)))


def Qform(Hm):
    """Q_0(M) = sum_{J not 0} (-1/2)^|J| Tr(Y_J^* M Y_J M) + sum_{J containing 0} (-1/2)^|J| ||Tr_J(C M)||^2"""
    val = 0
    for c, Y, pos, k, Jp in pieces:
        Mk = at_copy(Hm, pos, k); val += c * (Y.H * Mk * Y * Mk).trace()
        Z, _ = ptrace(Cs * at_copy(Hm, 0, 3), Jp + [0], 3); val += -c / 2 * sum(sp.expand(z * sp.conjugate(z)) for z in Z)
    return sp.nsimplify(sp.expand(val))


def Nform(Hm):
    Mk = at_copy(Hm, 0, 3); return sp.expand((Cs.H * Mk * Cs * Mk).trace())


Qm, Nm = sp.zeros(9, 9), sp.zeros(9, 9)
diagQ = [Qform(b) for b in basis]; diagN = [Nform(b) for b in basis]
for a in range(9):
    Qm[a, a], Nm[a, a] = diagQ[a], diagN[a]
    for b in range(a + 1, 9):
        Qm[a, b] = Qm[b, a] = sp.expand((Qform(basis[a] + basis[b]) - diagQ[a] - diagQ[b]) / 2)
        Nm[a, b] = Nm[b, a] = sp.expand((Nform(basis[a] + basis[b]) - diagN[a] - diagN[b]) / 2)
for lam in (sp.Rational(-1), sp.Rational(-1, 2)):
    minors = [(Qm - lam * Nm)[:k, :k].det() for k in range(1, 10)]
    good = all(m_ > 0 for m_ in minors); ok &= good
    print(f"n=3 exact: Q_0 - ({lam}) N_0 has all nine leading principal minors positive: {good} (smallest {min(minors)})")
# ---- n = 4, 5, 6: Gaussian-integer co-ranges ----
for n in (4, 5, 6):
    a = np.load(f"../data/sc_counterexample_corange_n{n}.npy", allow_pickle=False).astype(float); g = a[0] + 1j * a[1]
    u = tangent_range([np.eye(2)] * n, np.zeros(2), n)
    ok &= float_report(np.outer(u[0], g[0]) + np.outer(u[1], g[1]), n)
assert ok
print("PASS")
