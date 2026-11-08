"""Section 9: no degree-four certificate at three copies.  Degree-4 relaxation of TN_n (or of q_n >= 0) in the balanced
constant-diagonal gauge of Proposition 28.
Vectors v = (x1, x2, y1, y2) in H = (C^d)^n, C = x1 y1^* + x2 y2^*.  Data M_K[(a,b),(c,d)] = <Tr_K(a b^*), Tr_K(c d^*)>.
For a FIXED first-level Gram G1[a,c] = <a,c> (taken from an actual operator), minimise
    Delta = sum_K (-t)^|K| sum_ij M_K[(xi,yi),(xj,yj)] - m^{n-1}(N - t S^2)
over all M_K (K != empty, [n]) subject to:
  M_K PSD; adjoint symmetry M_K[(p,q),(r,s)] = M_K[(s,r),(q,p)]; crossing M_{K^c}[(a,b),(c,d)] = M_K[(d,b),(c,a)];
  sector Grams  <a(x)b|Pi_T|c(x)d> = 2^-n sum_L (-1)^{|L cap T|} M_{L^c}[(a,c),(d,b)]  PSD;
  component Grams  sum_{L >= K} (-1)^{|L\\K|} d^{-|L|} M_L  PSD.
M_empty[(a,b),(c,d)] = G1[a,c] G1[d,b] is fixed.
usage: python degree4_relaxation.py n t samples seed target [d]   (target = tn: the gap of (9.1); target = q: q_n itself)
At n = 2 the minimum is >= 0 (TN_2 holds, Proposition 28); at n = 3 it is about -0.08 S^2 at t = 1/2 for both targets,
although q_3 >= 0 is true: a proof at three copies needs information beyond these two-replica Gram data."""
import sys, itertools
import numpy as np
import cvxpy as cp

n, t, samples, seed, target = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
d = int(sys.argv[6]) if len(sys.argv) > 6 else 3
m = 1 - t
rng = np.random.default_rng(seed)
subsets = [frozenset(K) for r in range(n + 1) for K in itertools.combinations(range(n), r)]
full = frozenset(range(n))


def pm(fun):
    P = np.zeros((256, 256))
    for i in range(16):
        for j in range(16):
            i2, j2 = fun(i, j)
            P[j * 16 + i, j2 * 16 + i2] = 1
    return P


def dec(i):
    return divmod(i, 4)


P_cross = pm(lambda i, j: (4 * dec(j)[1] + dec(i)[1], 4 * dec(j)[0] + dec(i)[0]))   # E'[(a,b),(c,d)] = E[(d,b),(c,a)]
P_adj = pm(lambda i, j: (4 * dec(j)[1] + dec(j)[0], 4 * dec(i)[1] + dec(i)[0]))     # E'[(p,q),(r,s)] = E[(s,r),(q,p)]
P_F = pm(lambda i, j: (4 * dec(i)[0] + dec(j)[0], 4 * dec(j)[1] + dec(i)[1]))       # E'[(a,b),(c,d)] = E[(a,c),(d,b)]


def reidx(E, P):
    if isinstance(E, np.ndarray):
        return (P @ E.reshape(-1, order='F')).reshape((16, 16), order='F')
    return cp.reshape(P @ cp.vec(E, order='F'), (16, 16), order='F')


def gauge(C):
    U, s, Vh = np.linalg.svd(C)
    X0 = U[:, :2] * np.sqrt(s[:2]); Y0 = Vh.conj().T[:, :2] * np.sqrt(s[:2])
    M = Y0.conj().T @ X0
    M0 = M - np.trace(M) / 2 * np.eye(2)
    A = (M0 + M0.conj().T) / 2; B = (M0 - M0.conj().T) / 2j
    pa = np.array([A[0, 1].real, -A[0, 1].imag, A[0, 0].real]); pb = np.array([B[0, 1].real, -B[0, 1].imag, B[0, 0].real])
    r = np.cross(pa, pb)
    if np.linalg.norm(r) < 1e-12:
        r = np.cross(pa if np.linalg.norm(pa) > 1e-12 else np.array([1, 0, 0.]), np.array([0.3, 0.5, 0.7]))
    r /= np.linalg.norm(r)
    th, ph = np.arccos(np.clip(r[2], -1, 1)), np.arctan2(r[1], r[0])
    w = np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)])
    Q = np.column_stack([w, np.array([-np.conj(w[1]), np.conj(w[0])])])
    X, Y = X0 @ Q, Y0 @ Q
    return X, Y


def true_M(vs, K):
    D = vs[0].shape[0]
    out = np.zeros((16, 16), complex)
    trs = {}
    for a in range(4):
        for b in range(4):
            Z = np.outer(vs[a], vs[b].conj()).reshape([d] * (2 * n))
            Xr, rem = Z, list(range(n))
            for s_ in sorted(K, reverse=True):
                i = rem.index(s_); Xr = np.trace(Xr, axis1=i, axis2=len(rem) + i); rem.pop(i)
            trs[(a, b)] = np.asarray(Xr).reshape(-1)
    for i in range(16):
        for j in range(16):
            out[i, j] = np.vdot(trs[dec(i)], trs[dec(j)])
    return out


def run(C):
    X, Y = gauge(C)
    vs = [X[:, 0], X[:, 1], Y[:, 0], Y[:, 1]]
    G1 = np.array([[np.vdot(vs[a], vs[c]) for c in range(4)] for a in range(4)])
    Mfix = np.zeros((16, 16), complex)
    for i in range(16):
        for j in range(16):
            (a, b), (c, dd) = dec(i), dec(j)
            Mfix[i, j] = G1[a, c] * G1[dd, b]
    S = (G1[0, 0] + G1[1, 1]).real
    xy = [4 * 0 + 2, 4 * 1 + 3]                     # pairs (x1,y1), (x2,y2)
    N = sum(Mfix[i, j] for i in xy for j in xy).real
    Ms, cons = {}, []
    Ms[frozenset()] = Mfix; Ms[full] = reidx(Mfix, P_cross)
    lower = [K for K in subsets if K != frozenset() and K != full and sorted(K) < sorted(full - K)]
    for K in lower:
        V = cp.Variable((16, 16), hermitian=True)
        Ms[K] = V; Ms[full - K] = reidx(V, P_cross)
        cons += [V == reidx(V, P_adj)]
    for K in subsets:
        if not isinstance(Ms[K], np.ndarray):
            cons += [Ms[K] >> 0]
    for T in subsets:
        ST = sum(((-1) ** len(L & T)) * reidx(Ms[full - L], P_F) for L in subsets) / 2 ** n
        if not isinstance(ST, np.ndarray):
            cons += [(ST + ST.H) / 2 >> 0]
    for K in subsets:
        CK = sum(((-1) ** len(L - K)) * d ** (-len(L)) * Ms[L] for L in subsets if K <= L)
        if not isinstance(CK, np.ndarray):
            cons += [(CK + CK.H) / 2 >> 0]
    qn = sum(((-t) ** len(K)) * sum(Ms[K][i, j] for i in xy for j in xy) for K in subsets)
    obj = cp.real(qn) - (m ** (n - 1) * (N - t * S * S) if target == 'tn' else 0)
    prob = cp.Problem(cp.Minimize(obj), cons)
    prob.solve(solver='CLARABEL')
    # true value at the actual vectors
    true = sum(((-t) ** len(K)) * sum(true_M(vs, K)[i, j] for i in xy for j in xy) for K in subsets).real - (m ** (n - 1) * (N - t * S * S) if target == 'tn' else 0)
    return prob.value / S ** 2, true / S ** 2, prob.status


for k in range(samples):
    D = d ** n
    kind = k % 3
    if kind == 0:     # generic complex rank two
        C = (rng.normal(size=(D, 2)) + 1j * rng.normal(size=(D, 2))) @ (rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D)))
    elif kind == 1:   # near the product family X (x) |0><0|^{n-1}
        Xp = rng.normal(size=(d, 2)); Xp = Xp @ Xp.T
        e0 = np.zeros(d ** (n - 1)); e0[0] = 1
        C = np.kron(Xp, np.outer(e0, e0)) + 0.05 * (rng.normal(size=(D, 2)) @ rng.normal(size=(2, D)))
        U, s, Vh = np.linalg.svd(C); C = U[:, :2] @ np.diag(s[:2]) @ Vh[:2]
    else:             # near the signed GHZ pair
        e0 = np.zeros(D); e0[0] = 1; e1 = np.zeros(D); e1[int('1' * n, d)] = 1
        C = np.outer(e0, e0) - (-1) ** n * np.outer(e1, e1) + 0.05 * (rng.normal(size=(D, 2)) @ rng.normal(size=(2, D)))
        U, s, Vh = np.linalg.svd(C); C = U[:, :2] @ np.diag(s[:2]) @ Vh[:2]
    val, true, st = run(C)
    print(f"n={n} t={t} target={target} sample {k} ({['generic', 'near product', 'near GHZ'][kind]}): relaxation min/S^2 = {val:+.6f}   true value/S^2 = {true:+.6f}   [{st}]", flush=True)
