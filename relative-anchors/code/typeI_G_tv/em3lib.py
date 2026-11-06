"""The third LDL^T pivot of (B_inv) with both level sizes symbolic, by the determinant route (NOTES.md, Section 10).

The matrix.  For s_t = s (mode T11) or s_t = 1 (mode T22), A is the 3 x 3 matrix of (B_inv) exactly as in
prove_gersh.run_inv (y-coordinates, H = 1/2):
    A = G_y - H (1 - theta_t) diag(kappa),   G_y = N^T G_z N,   kappa = (1, s_t sigma/(S (sigma + S)), s_t S1 S2/S),
computed from the validated closed form tl_inv_gen_T11/T22 in its own variables (u1, u2, sig, s1), and then composed,
entry by entry, into the orthant variables (x1, x2, r, s) by u_i = 2 + x_i, sigma = s_t + r, s1 = s (an invertible
affine change of variables, so this is the same rational function as the late composition of prove_gersh.py).
A > 0 at a point iff its three leading principal minors are positive; the first two pivots are certified by
prove_gersh.py (gen, B_s12 and B_112), so the third pivot is positive iff det A > 0.

The bookkeeping (as in work_g2/g2lib.py and work_r2b/f2det.py):
  * every irreducible factor f of an entry denominator of A, its multiplicity m_f(a, b) in den A_ab, and
    e_f = max over permutations p of sum_a m_f(a, p(a)) (exact dynamic programming over subsets);
  * D* = prod f^e_f; then T = D* det A is a polynomial (each Leibniz term D* prod_a A_{a,p(a)} is one);
  * deg_x1 T <= B = deg_x1 D* + max_p sum_a (deg_x1 num A_{a,p(a)} - deg_x1 den A_{a,p(a)}).
Slices: A with x1 = c fixed, exact determinant over Q(x2, r, s)."""
import os, sys, time, pickle, gc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
import prove_gersh as PG
from flintfield import Field, F as FE

Q = flint.fmpq
VN = PG.VN                      # ('x1', 'x2', 'r', 's')
KX = PG.KX
ONE = PG.ONE
MODE = {'s': 'T11', '1': 'T22'}
SLICE_VARS = ('x2', 'r', 's')
K3 = Field(list(SLICE_VARS))
ONE3 = K3.ctx.from_dict({(0, 0, 0): 1})


def build_A(ka, shifts=(2, 2), H=Q(1, 2), log=print, orig=False):
    """A (3 x 3, list of lists of flintfield fractions over Q(x1, x2, r, s)) for s_t = ka ('s' or '1'), with
    u1 = shifts[0] + x1, u2 = shifts[1] + x2, sigma = s_t + r, s1 = s.  The algebra is that of prove_gersh.run_inv
    (non-bz branch), line by line."""
    t0 = time.time()
    names, K, GK, gen = PG.orig('tl_inv_gen', MODE[ka])
    one = K(1)
    wK = {'s': gen['s1'], '1': one}
    a, sig = wK[ka], gen['sig']
    S1, S2 = PG.level_sums(K, gen, None, None)
    S = S1 + S2
    Gz = [[GK[(min(i, j), max(i, j))] for j in range(3)] for i in range(3)]
    del GK
    om_th = one * H * (sig - a) / sig
    Nm = [[one, K(0), K(0)], [K(0), one / S, S2 / S], [K(0), one / S, -S1 / S]]
    T = [[sum((Gz[k][l] * Nm[l][j] for l in range(3) if Nm[l][j] != 0), K(0)) for j in range(3)] for k in range(3)]
    Gy = [[sum((Nm[k][i] * T[k][j] for k in range(3) if Nm[k][i] != 0), K(0)) for j in range(3)] for i in range(3)]
    kap = [one, a * sig / (S * (sig + S)), a * S1 * S2 / S]
    A = [[Gy[i][j] - (om_th * kap[i] if i == j else K(0)) for j in range(3)] for i in range(3)]
    del T, Gy, Gz
    gc.collect()
    sym = all(A[i][j] == A[j][i] for i in range(3) for j in range(i + 1, 3))
    assert sym
    x1, x2, r, s = (PG.GX[v] for v in VN)
    wts_x = {'s': s, '1': ONE}
    mp = dict(u1=shifts[0] * ONE + x1, u2=shifts[1] * ONE + x2, sig=wts_x[ka] + r, s1=s)
    Ax = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(i, 3):
            Ax[i][j] = PG.late(A[i][j], names, mp)
            Ax[j][i] = Ax[i][j]
    log('  A built for s_t = %s (u_i = %d + x_i): symmetric %s; entry numerator/denominator terms %s [%.1fs]'
        % (ka, shifts[0], sym, [(len(Ax[i][j].n), len(Ax[i][j].d)) for i in range(3) for j in range(i, 3)],
           time.time() - t0), flush=True)
    if orig:
        return Ax, (names, A)
    return Ax


def maxperm(W):
    """max over permutations p of sum_a W[a][p(a)] (W[a][b] = None: forbidden); exact DP over subsets"""
    n = len(W)
    best = {0: 0}
    for a in range(n):
        nb = {}
        for mask, val in best.items():
            for b in range(n):
                if mask >> b & 1 or W[a][b] is None:
                    continue
                m2 = mask | (1 << b)
                v = val + W[a][b]
                if m2 not in nb or v > nb[m2]:
                    nb[m2] = v
        best = nb
    return best.get((1 << n) - 1)


def factor_data(A):
    """(factors, mult, estar): the irreducible factors of the entry denominators (in a fixed order), the
    multiplicities mult[(k, a, b)], and the exponents e_f of D*"""
    keys, facs, mult = [], {}, {}
    n = len(A)
    for a in range(n):
        for b in range(n):
            if A[a][b].n.is_zero():
                continue
            c, fl = A[a][b].d.factor()
            assert c > 0
            for f, m in fl:
                k = str(f)
                if k not in facs:
                    facs[k] = f
                    keys.append(k)
                mult[(keys.index(k), a, b)] = m
    factors = [facs[k] for k in keys]
    estar = []
    for kf in range(len(factors)):
        W = [[(mult.get((kf, a, b), 0) if not A[a][b].n.is_zero() else None) for b in range(n)] for a in range(n)]
        estar.append(maxperm(W))
    return factors, mult, estar


def degree_bound(A, factors, estar, iv):
    """B_v = deg_v D* + max_p sum_a (deg_v num - deg_v den), v = VN[iv]: an upper bound for deg_v (D* det A)"""
    n = len(A)
    W = [[(A[a][b].n.degrees()[iv] - A[a][b].d.degrees()[iv]) if not A[a][b].n.is_zero() else None
          for b in range(n)] for a in range(n)]
    dD = sum(e * f.degrees()[iv] for f, e in zip(factors, estar))
    return dD + maxperm(W), dD


def spec(P, c):
    """the polynomial P(x1, x2, r, s) at x1 = c, in the context of (x2, r, s)"""
    g = K3.ctx.gens()
    return P.compose(K3.ctx.from_dict({(0, 0, 0): Q(c)}), g[0], g[1], g[2], ctx=K3.ctx)


def det3(M):
    """exact determinant of a 3 x 3 symmetric matrix of fractions: det = M00 det(Schur complement of M00)
    (M00 is nonzero on the region)"""
    p1 = M[0][0]
    f1 = M[0][1] / p1
    f2 = M[0][2] / p1
    b11 = M[1][1] - f1 * M[0][1]
    b12 = M[1][2] - f1 * M[0][2]
    b22 = M[2][2] - f2 * M[0][2]
    return p1 * (b11 * b22 - b12 * b12)


def slice_det(A, c):
    """det A at x1 = c, a reduced fraction over Q(x2, r, s), computed in the orthant variables"""
    M = [[FE(K3, spec(A[i][j].n, c), spec(A[i][j].d, c)) for j in range(3)] for i in range(3)]
    return det3(M)


class OrigSlicer:
    """det A at x1 = c computed in the variables of the closed form: the entries of A over Q(u1, u2, sig, s1) are
    specialised at u1 = shift1 + c, the determinant is taken over Q(u2, sig, s1), and only the result is composed
    into (x2, r, s) by u2 = shift2 + x2, sig = s_t + r, s1 = s (an invertible affine change of variables, so the
    result is the same reduced fraction as slice_det)."""

    def __init__(self, ka, names, Aorig, shifts=(2, 2)):
        assert tuple(names) == ('u1', 'u2', 'sig', 's1')
        self.ka, self.A, self.shifts = ka, Aorig, shifts
        self.K = Field(['u2', 'sig', 's1'])
        g = K3.ctx.gens()                                 # x2, r, s
        x2, r, s = g
        cst = lambda v: K3.ctx.from_dict({(0, 0, 0): Q(v)})
        self.back = [cst(shifts[1]) + x2, (s if ka == 's' else cst(1)) + r, s]

    def spec_u1(self, P, N):
        g = self.K.ctx.gens()
        return P.compose(self.K.ctx.from_dict({(0, 0, 0): Q(N)}), g[0], g[1], g[2], ctx=self.K.ctx)

    def det(self, c):
        N = self.shifts[0] + c
        M = [[None] * 3 for _ in range(3)]
        for i in range(3):
            for j in range(i, 3):
                M[i][j] = FE(self.K, self.spec_u1(self.A[i][j].n, N), self.spec_u1(self.A[i][j].d, N))
                M[j][i] = M[i][j]
        d = det3(M)
        return FE(K3, d.n.compose(*self.back, ctx=K3.ctx), d.d.compose(*self.back, ctx=K3.ctx))


def multiplicity(f, P):
    """the largest e with f^e | P (f non-constant)"""
    e = 0
    while True:
        q, r = divmod(P, f)
        if r != 0:
            return e
        P = q
        e += 1


def one_signed_pos(f, free):
    """True if every coefficient of f is > 0 and f has a monomial free of the variables `free` (so f > 0 for
    x1, x2, r >= 0 and s > 0 when free = (x1, x2, r))"""
    idx = [VN.index(v) for v in free]
    return all(c > 0 for c in f.coeffs()) and any(all(m[i] == 0 for i in idx) for m in f.monoms())
