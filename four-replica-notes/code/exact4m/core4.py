"""Exact data for the two-letter relaxation at n = 4 (d = 3), built for modular arithmetic.

Conventions follow ../relaxations/exact_coreN.py (see ../../NOTES.md, Section 3):
  * per-site algebras A (S_4 words), B1 (B_{3,1}(3), realignment target), A3 (S_3 words) with exact rational irreps
    (exact_site.ExactAlg): block l has dimension f_l, multiplicity m_l, Gram matrix G_l and rho_l(sigma);
  * variables: for each family c and sorted site pattern lam, B_lam = V Y V^T (V an exact basis of range(Pi^T), Pi the
    average over the letter stabiliser of (x)_s rho_{lam_s}(sigma)), Y symmetric and invariant under the site
    stabiliser of lam; X_hat_lam = G_lam^{-1} B_lam;  X_c >= 0  <=>  Y >= 0 for every (c, lam);
  * here the basis Y_k of the invariant symmetric matrices is the reduced row echelon basis, so the coordinate of an
    invariant Y on Y_k is the entry of Y at the pivot position of Y_k.
Integer scaling: for a group (c, lam) with K parameters, Xt = dX * X_hat (K, D, D) is an integer tensor (Python ints).
"""
import sys, os, time, itertools, pickle, math, numpy as np, flint
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'relaxations'))
from exact_core import site_algebras, FAMS, letter_group, locate, exact_coeffs, SiteMap
from exact_site import nullspace_q
from alg4 import S4, PERM, ptrace_last, inv
from rep4 import realign1

NS = 4                                   # number of sites
fq = flint.fmpq

# ----------------------------------------------------------------------------- small exact helpers
def fm_to_frac(M):
    """fmpq_mat -> numpy object array of (num, den) python Fractions"""
    from fractions import Fraction
    return np.array([[Fraction(int(M[i, j].p), int(M[i, j].q)) for j in range(M.ncols())] for i in range(M.nrows())],
                    dtype=object).reshape(M.nrows(), M.ncols())

def fm_from_frac(A):
    A = np.asarray(A, dtype=object)
    from fractions import Fraction
    vals = []
    for x in A.ravel():
        x = Fraction(x); vals.append(fq(x.numerator, x.denominator))
    return flint.fmpq_mat(A.shape[0], A.shape[1], vals)

def fm_kron(mats):
    """Kronecker product of fmpq_mats (python loops; sizes here are <= 81)"""
    M = flint.fmpq_mat(1, 1, [1])
    for X in mats:
        ra, ca, rb, cb = M.nrows(), M.ncols(), X.nrows(), X.ncols()
        R = flint.fmpq_mat(ra * rb, ca * cb)
        for i in range(ra):
            for j in range(ca):
                a = M[i, j]
                if a == 0: continue
                for k in range(rb):
                    for l in range(cb):
                        b = X[k, l]
                        if b != 0: R[i * rb + k, j * cb + l] = a * b
        M = R
    return M

def fm_den(M):
    """lcm of the denominators of an fmpq_mat"""
    d = 1
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            q = int(M[i, j].q)
            if q != 1: d = d * q // math.gcd(d, q)
    return d

def fm_to_int(M, d):
    """integer numpy object array d * M (d a common denominator)"""
    out = np.empty((M.nrows(), M.ncols()), dtype=object)
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            q = M[i, j]; v = int(q.p) * d
            assert v % int(q.q) == 0
            out[i, j] = v // int(q.q)
    return out

def perm_matrix_fm(shape, p):
    """P with P e_(i_0..i_{n-1}) = e_(i'), i'_q = i_{p[q]}  (site permutation of tensor factors)"""
    D = int(np.prod(shape)); n = len(shape)
    P = flint.fmpq_mat(D, D)
    for idx in itertools.product(*[range(s) for s in shape]):
        a = np.ravel_multi_index(idx, shape)
        new = tuple(idx[p[q]] for q in range(n)); nshape = tuple(shape[p[q]] for q in range(n))
        b = np.ravel_multi_index(new, nshape)
        P[b, a] = 1
    return P

def sorted_patterns(nb, n=NS):
    return list(itertools.combinations_with_replacement(range(nb), n))

def orbit_perms(lam):
    """distinct site permutations of a sorted pattern: list of (lam', perm), lam'[q] = lam[perm[q]]"""
    n = len(lam); seen = {}
    for perm in itertools.permutations(range(n)):
        lp = tuple(lam[perm[q]] for q in range(n))
        if lp not in seen: seen[lp] = perm
    return [(lp, p) for lp, p in seen.items()]

def stabilizer(lam):
    n = len(lam)
    return [p for p in itertools.permutations(range(n)) if all(lam[p[q]] == lam[q] for q in range(n))]


# ----------------------------------------------------------------------------- variable space
class VarSpace4:
    """exact variable space at n = 4; groups[(c, lam)] = dict(V, Ys (list of fmpq_mat), piv (pivot (a, b) per Y_k),
    idx, t, D, shape, dX, Xt (object int array (K, D, D)), G (fmpq_mat), Gi (fmpq_mat))"""
    def __init__(self, A, log=print, with_tensors=True):
        t0 = time.time(); self.A = A; nb = len(A.blocks); n = NS
        self.f = [b['f'] for b in A.blocks]; self.m = [b['m'] for b in A.blocks]
        self.N = 0; self.groups = {}; self.order = []
        for c in FAMS:
            Gc = letter_group(c)
            for lam in sorted_patterns(nb):
                shape = tuple(self.f[l] for l in lam); D = int(np.prod(shape))
                Pi = flint.fmpq_mat(D, D)
                for s in Gc: Pi = Pi + fm_kron([A.blocks[l]['rho'][s] for l in lam])
                Pi = Pi * fq(1, len(Gc))
                R, t = Pi.rref()
                if t == 0: continue
                V = flint.fmpq_mat(D, t, [R[k, j] for j in range(D) for k in range(t)])     # columns = RREF rows
                stab = [p for p in stabilizer(lam) if list(p) != list(range(n))]
                VtV = V.transpose() * V; VtVi = VtV.inv()
                Ks = []
                for p in stab:
                    Pq = perm_matrix_fm(shape, p)
                    K = VtVi * (V.transpose() * (Pq * V))
                    assert Pq * V == V * K, "site permutation does not preserve range(Pi^T)"
                    Ks.append(K)
                # invariant symmetric matrices: Reynolds images of the symmetric units, RREF basis (as t*t vectors)
                pairs = [(a, b_) for a in range(t) for b_ in range(a, t)]
                rows = []
                for (a, b_) in pairs:
                    E = flint.fmpq_mat(t, t); E[a, b_] = 1; E[b_, a] = 1
                    Y = E
                    for K in Ks: Y = Y + K * E * K.transpose()
                    rows.append([Y[i, j] for i in range(t) for j in range(t)])
                Mrow = flint.fmpq_mat(len(rows), t * t, [v for r in rows for v in r])
                RR, kY = Mrow.rref()
                Ys = []; piv = []
                for k in range(kY):
                    Y = flint.fmpq_mat(t, t, [RR[k, e] for e in range(t * t)])
                    # pivot = first nonzero entry of the row
                    e0 = next(e for e in range(t * t) if RR[k, e] != 0)
                    assert RR[k, e0] == 1
                    Ys.append(Y); piv.append(divmod(e0, t))
                    assert Y == Y.transpose()
                idx = list(range(self.N, self.N + kY)); self.N += kY
                Gl = fm_kron([A.blocks[l]['G'] for l in lam]); Gil = fm_kron([A.blocks[l]['Gi'] for l in lam])
                g = dict(c=c, lam=lam, V=V, Ys=Ys, piv=piv, idx=idx, t=t, D=D, shape=shape, G=Gl, Gi=Gil, K=kY,
                         orbit=orbit_perms(lam), mlam=int(np.prod([self.m[l] for l in lam])))
                if with_tensors:
                    Xh = [Gil * (V * Y * V.transpose()) for Y in Ys]
                    dX = 1
                    for M in Xh:
                        dd = fm_den(M); dX = dX * dd // math.gcd(dX, dd)
                    g['dX'] = dX
                    g['Xt'] = np.array([fm_to_int(M, dX) for M in Xh], dtype=object)      # (K, D, D)
                self.groups[(c, lam)] = g; self.order.append((c, lam))
        log(f"exact variable space n=4: {self.N} parameters, {len(self.groups)} groups; {time.time() - t0:.1f}s")

    def pair_layout(self, X, shape):
        """(K, D, D) -> (K, f0*f0, f1*f1, f2*f2, f3*f3) with pair index (i_s, j_s)"""
        K = X.shape[0]; n = len(shape)
        Y = X.reshape((K,) + shape + shape)
        Y = Y.transpose([0] + [a for q in range(n) for a in (1 + q, 1 + n + q)])
        return Y.reshape((K,) + tuple(s * s for s in shape))
