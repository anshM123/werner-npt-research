"""Checks for Theorem 7.19 (tensor products of plane maps with different parameters, every n).

Q_a(Y) = a Tr(Y) I_2 - (Y + Y^T)/4 on a real plane, Q_a = (2a - 1/2) Ad_Y o Y_t with t = 1/(4a - 1), where Y_t is the
Pauli map with eigenvalues (1, t, 0, t) on (I, X, Y, Z).  Write tau = 2t = 1/(2a - 1/2).  The theorem says that
Q_{a_1} (x) ... (x) Q_{a_n} is positive, and then decomposable, exactly when prod_{s in S} (2 a_s - 1/2) >= 1/2 for every
nonempty S, i.e. prod_{s in S} tau_s <= 2.

1. exact (fractions): for rational tau_s >= 1 with prod tau_s = 2, the Bell-diagonal p = (F + F~)/2 - R, with
   F = prod (1 + tau_s, 1, 1, 1 - tau_s), F~ = prod (1 + tau_s, 1, 1, tau_s - 1), R = prod (2/tau_s, 1, 1, 0),
   satisfies p + gamma^(x)n p = F and p >= 0 (gamma = partial transpose on Bell-diagonal operators).
2. explicit matrices (floating point, n <= 3): P >= 0 and P + P^Gamma = (x)_s C(Y_{t_s}), partial transpose on the matrix.
3. control: with prod tau_s = 2.1 the same p does not reproduce F.
4. witnesses: if the condition fails on S, the plane GHZ vector on S (product vectors elsewhere) gives
   <y| (x)_s Q_{a_s}(x x^T) |y> = 2^-|S| [prod_S (2 a_s - 1/2) - 1/2] prod_{s not in S} a_s < 0 (explicit vectors).
5. consistency (floating point, n = 3): alternating minimisation of the least output eigenvalue over real plane inputs at
   random parameter vectors agrees in sign with the condition.
Ends with PASS or FAIL."""
import itertools
from fractions import Fraction as Fr
import numpy as np

ok = True
gam = [[Fr(1, 2) * x for x in row] for row in [[1, 1, 1, -1], [1, 1, -1, 1], [1, -1, 1, 1], [-1, 1, 1, 1]]]

# 1. exact identity and positivity for rational parameter vectors with product 2
cases = [[Fr(2)], [Fr(4, 3), Fr(3, 2)], [Fr(5, 4), Fr(8, 5)], [Fr(1), Fr(2)], [Fr(6, 5), Fr(5, 4), Fr(4, 3)],
         [Fr(1), Fr(8, 7), Fr(7, 4)], [Fr(8, 7), Fr(7, 6), Fr(6, 5), Fr(5, 4)], [Fr(10, 9), Fr(9, 8), Fr(4, 3), Fr(6, 5)]]
for taus in cases:
    n = len(taus)
    assert np.prod([float(t) for t in taus]) == 2.0 or abs(np.prod([float(t) for t in taus]) - 2) < 1e-12
    fv = [[1 + t, Fr(1), Fr(1), 1 - t] for t in taus]
    ftv = [[1 + t, Fr(1), Fr(1), t - 1] for t in taus]
    rv = [[2 / t, Fr(1), Fr(1), Fr(0)] for t in taus]
    def pr(vecs, J):
        out = Fr(1)
        for s, j in enumerate(J): out *= vecs[s][j]
        return out
    strings = list(itertools.product(range(4), repeat=n))
    p = {J: Fr(1, 2) * (pr(fv, J) + pr(ftv, J)) - pr(rv, J) for J in strings}
    good = True
    for J in strings:
        gp = Fr(0)
        for K in strings:
            c = Fr(1)
            for s in range(n): c *= gam[J[s]][K[s]]
            if c: gp += c * p[K]
        good &= (p[J] + gp == pr(fv, J))
    pmin = min(p.values())
    c_ok = good and pmin >= 0
    print(f"1. tau = {[str(t) for t in taus]}: p + gamma p = F exactly: {good}; min p = {pmin} ({float(pmin):.4f})")
    ok &= c_ok

# 2. explicit matrices
bell = [np.array([1, 0, 0, 1]) / np.sqrt(2), np.array([1, 0, 0, -1]) / np.sqrt(2),
        np.array([0, 1, 1, 0]) / np.sqrt(2), np.array([0, 1, -1, 0]) / np.sqrt(2)]
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Yp = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1])
def ptrans_all(M, n):
    T = M.reshape([2] * (4 * n)); perm = list(range(4 * n))
    for s in range(n):
        rb, cb = 2 * s + 1, 2 * n + 2 * s + 1
        perm[rb], perm[cb] = perm[cb], perm[rb]
    return T.transpose(perm).reshape(4**n, 4**n)
for taus in [[2.0], [4 / 3, 1.5], [1.25, 1.6], [1.2, 1.25, 4 / 3], [1.0, 8 / 7, 1.75]]:
    n = len(taus)
    Cn = np.ones((1, 1))
    for tau in taus:
        t = tau / 2
        C1 = 0.25 * sum(l * np.kron(P, P.conj()) for l, P in zip([1, t, 0, t], (I2, X, Yp, Z)))
        Cn = np.kron(Cn, C1)
    Pn = np.zeros((4**n, 4**n)); mn = np.inf
    for J in itertools.product(range(4), repeat=n):
        f = np.prod([[1 + tt, 1, 1, 1 - tt][j] for tt, j in zip(taus, J)])
        ft = np.prod([[1 + tt, 1, 1, tt - 1][j] for tt, j in zip(taus, J)])
        r = np.prod([[2 / tt, 1, 1, 0][j] for tt, j in zip(taus, J)])
        val = 0.5 * (f + ft) - r; mn = min(mn, val)
        v = bell[J[0]]
        for j in J[1:]: v = np.kron(v, bell[j])
        Pn += val / 4**n * np.outer(v, v)
    err = np.abs(Pn + ptrans_all(Pn, n) - Cn).max()
    m_ok = mn > -1e-12 and err < 1e-12
    print(f"2. tau = {np.round(taus, 4)}: min eigenvalue of 4^n P = {mn:.2e}; max |P + P^Gamma - C| = {err:.1e}")
    ok &= m_ok

# 3. control: product 2.1
taus = [Fr(3, 2), Fr(7, 5)]          # product 21/10
n = 2
fv = [[1 + t, Fr(1), Fr(1), 1 - t] for t in taus]; ftv = [[1 + t, Fr(1), Fr(1), t - 1] for t in taus]
rv = [[2 / t, Fr(1), Fr(1), Fr(0)] for t in taus]
def pr(vecs, J):
    out = Fr(1)
    for s, j in enumerate(J): out *= vecs[s][j]
    return out
strings = list(itertools.product(range(4), repeat=n))
p = {J: Fr(1, 2) * (pr(fv, J) + pr(ftv, J)) - pr(rv, J) for J in strings}
res = max(abs(p[J] + sum(gam[J[0]][K[0]] * gam[J[1]][K[1]] * p[K] for K in strings) - pr(fv, J)) for J in strings)
c_ok = res > 0
print(f"3. control, prod tau = 21/10: identity fails (max residual {float(res):.4f}): {c_ok}")
ok &= c_ok

# 4. witnesses with explicit vectors
def apply_Q(a_list, Xm):
    n = len(a_list); T = Xm.reshape([2] * (2 * n)).astype(complex)
    for s, a in enumerate(a_list):
        T = np.moveaxis(T, [s, n + s], [0, 1]); sh = T.shape
        Mi = T.reshape(2, 2, -1); out = np.zeros_like(Mi); tr = np.einsum('iir->r', Mi)
        for k in range(2): out[k, k] += a * tr
        out -= 0.25 * (Mi + Mi.transpose(1, 0, 2))
        T = out.reshape(sh); T = np.moveaxis(T, [0, 1], [s, n + s])
    return T.reshape(2**n, 2**n)
pv = np.array([1, 1j]) / np.sqrt(2); e = np.array([1.0, 0.0]); eperp = np.array([0.0, 1.0])
Ymat = np.array([[0, -1.0], [1.0, 0]])
for a_list, S in [([0.6, 0.6, 0.9], [0, 1]), ([0.6, 0.64, 0.66], [0, 1, 2]), ([0.55, 0.66, 0.8, 1.0], [0, 1]),
                  ([0.64, 0.65, 0.66, 0.67], [0, 1, 2, 3]), ([0.45, 0.9, 0.9], [0])]:
    n = len(a_list)
    g1 = np.ones(1, dtype=complex); g2 = np.ones(1, dtype=complex)
    for s in range(n):
        if s in S: g1 = np.kron(g1, pv); g2 = np.kron(g2, pv.conj())
        else: g1 = np.kron(g1, e); g2 = np.kron(g2, e)
    x = ((g1 + g2) / np.sqrt(2)).real
    if len(S) == 1:
        x = np.ones(1)
        for s in range(n): x = np.kron(x, e)
        y = np.ones(1)
        for s in range(n): y = np.kron(y, e if s in S else eperp)
    else:
        y = x.copy()
        if len(S) % 2 == 0:   # Y on the first copy of S
            ops = [Ymat if s == S[0] else np.eye(2) for s in range(n)]
            Op = ops[0]
            for o in ops[1:]: Op = np.kron(Op, o)
            y = Op @ y
        # outputs outside S: e_perp instead of e
        ops = [np.eye(2) if s in S else np.array([[0, -1.0], [1.0, 0]]) for s in range(n)]
        Op = ops[0]
        for o in ops[1:]: Op = np.kron(Op, o)
        y = Op @ y
    val = np.real(y @ apply_Q(a_list, np.outer(x, x)) @ y)
    pred = 2.0**(-len(S)) * (np.prod([2 * a_list[s] - 0.5 for s in S]) - 0.5) * np.prod([a_list[s] for s in range(n) if s not in S])
    w_ok = abs(val - pred) < 1e-12 and val < 0
    print(f"4. a = {a_list}, S = {S}: value {val:+.6f}, formula {pred:+.6f}: {w_ok}")
    ok &= w_ok

# 5. consistency of signs at random parameter vectors, n = 3
rng = np.random.default_rng(7)
def minval(a_list, starts=8, iters=60):
    n = len(a_list); best = np.inf
    for st in range(starts):
        x = rng.standard_normal(2**n); x /= np.linalg.norm(x)
        for it in range(iters):
            W = apply_Q(a_list, np.outer(x, x)).real; w, V = np.linalg.eigh((W + W.T) / 2); y = V[:, 0]
            W2 = apply_Q(a_list, np.outer(y, y)).real; w2, V2 = np.linalg.eigh((W2 + W2.T) / 2); x = V2[:, 0]
        best = min(best, w2[0])
    # also the witnesses
    for r in range(1, n + 1):
        for S in itertools.combinations(range(n), r):
            best = min(best, 2.0**(-r) * (np.prod([2 * a_list[s] - 0.5 for s in S]) - 0.5) * np.prod([a_list[s] for s in range(n) if s not in S]))
    return best
agree = 0; tot = 0
for _ in range(40):
    a_list = list(rng.uniform(0.5, 0.85, 3))
    cond = min(np.prod([2 * a_list[s] - 0.5 for s in S]) for r in range(1, 4) for S in itertools.combinations(range(3), r)) >= 0.5
    mv = minval(a_list)
    tot += 1; agree += ((mv >= -1e-10) == cond)
s_ok = agree == tot
print(f"5. random parameter vectors (n = 3): sign of the least value agrees with the condition in {agree}/{tot}")
ok &= s_ok
print("PASS" if ok else "FAIL")
