"""General Pauli maps (eigenvalues 1, l1, l2, l3 on I, X, Y, Z).  Compare, for n copies:
 (a) FM's GHZ necessary condition (their Proposition 7, eq. (38) of arXiv v3);
 (b) the canonical certificate: P = E_n(C) - (Gamma|Gamma C|)^n / 2,  Q = E_n(Gamma C) - (Gamma|C|)^n / 2,
     E_n(M) = (M^n + |M|^n)/2, which satisfies P + Gamma Q = C^n identically; decomposable if P, Q >= 0 (pointwise on
     Bell strings, everything being Bell-diagonal);
 (c) decomposability by the exact LP (reduced to count classes is not valid without symmetry; use all 4^n strings).
Bell labels Phi+, Phi-, Psi+, Psi-; eigenvalue of sum lam_P P(x)conj(P)/4 on them: chars below."""
import numpy as np, itertools, sys
from scipy.optimize import linprog
chi = np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, -1, -1, 1], [1, 1, -1, -1]], float)   # rows I,X,Y,Z; cols Bell
gam = 0.5 * np.array([[1, 1, 1, -1], [1, 1, -1, 1], [1, -1, 1, 1], [-1, 1, 1, 1]], float)
def bell(l):  # Bell eigenvalues x4 of the Choi matrix
    lam = np.array([1, l[0], l[1], l[2]])
    return lam @ chi
def fm_ok(l, n):
    for (i, j, k) in itertools.permutations(range(3)):
        for p_ in range(n + 1):
            q_ = n - p_
            lhs = abs((1 + l[i])**p_ * (1 - l[i])**q_ + (1 - l[i])**p_ * (1 + l[i])**q_)
            for sgn in (1, -1):
                rhs = abs((l[j] + l[k])**p_ * (l[j] - l[k])**q_ + sgn**n * (l[j] - l[k])**p_ * (l[j] + l[k])**q_)
                if lhs < rhs - 1e-12: return False
    return True
def canon_ok(l, n):
    f = bell(l); fp = gam @ f                       # Gamma C
    af, afp = np.abs(f), np.abs(fp)
    g1, g2 = gam @ afp, gam @ af                   # Gamma|Gamma C|, Gamma|C|
    for J in itertools.product(range(4), repeat=n):
        J = list(J)
        P = 0.5 * (np.prod(f[J]) + np.prod(af[J])) - 0.5 * np.prod(g1[J])
        Q = 0.5 * (np.prod(fp[J]) + np.prod(afp[J])) - 0.5 * np.prod(g2[J])
        if P < -1e-12 or Q < -1e-12: return False
    return True
def lp_ok(l, n):
    f = bell(l); strings = list(itertools.product(range(4), repeat=n)); N = len(strings)
    F = np.array([np.prod(f[list(J)]) for J in strings])
    Gn = np.ones((N, N))
    for a_, J in enumerate(strings):
        for b_, K in enumerate(strings):
            Gn[a_, b_] = np.prod([gam[J[s], K[s]] for s in range(n)])
    # find q >= 0 with F - Gn q >= 0
    res = linprog(np.zeros(N), A_ub=Gn, b_ub=F, bounds=[(0, None)] * N, method='highs')
    return res.status == 0
def positive_n(l, n, trials=0):
    return None
if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    n = int(sys.argv[1])
    cnt = {}
    for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 300):
        l = rng.uniform(-1, 1, 3)
        f = bell(l); fp = gam @ f
        if f.min() >= 0 or fp.min() >= 0:      # CP or coCP: trivial
            continue
        a_, b_, c_ = fm_ok(l, n), canon_ok(l, n), lp_ok(l, n)
        key = (a_, b_, c_); cnt[key] = cnt.get(key, 0) + 1
        if a_ and not b_:
            if cnt[key] <= 3: print('FM ok but canonical fails:', np.round(l, 4), ' LP decomposable:', c_)
    print(f'n = {n}: counts (FM necessary, canonical certificate, LP decomposable): {cnt}')
