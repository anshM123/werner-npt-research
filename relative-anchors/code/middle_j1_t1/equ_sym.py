"""Equal U-weights, every k: the P-first reduction Phi of the |T| = 1 block of (R_inf)_1, symbolically in the number
u = k - 1 of copies in U.

Setting.  T = {t}, U = [k] minus t, u = |U| >= 2 (k >= 3).  Unnormalised squared weights: s_t = p, s_h = q (h in U),
L = p + u q; normalised x_t = sigma = p/L, x_h = s = q/L.  Every closed form of tlib.py (Lemma T1.3) is used with these
x's (the unrescaled block is homogeneous of degree 0 in the weights).

Rescaled frame.  D = diag(mu) with mu = w_c on R_c, P_c, Q_c.  For an operator O between coordinate spaces R^[k] or
R^U put O~ = D O D (entries mu_i O_ij mu_j, rational in p, q).  Then (O1 O2)~ = O1~ D^{-2} O2~,
(O^{-1})~ = D^2 (O~)^{-1} D^2 and (O^T)~ = (O~)^T.

Symmetry.  The permutations of U act on everything, and D commutes with them.  An invariant operator is stored by its
matrix in the bases (e_t, 1_U) of R^[k] and (1_U) of R^U, and by its scalar on the standard part {f on U : sum f = 0}.
For a k x k invariant operator with entries O[t,t] = al, O[t,h] = be, O[h,t] = bp, O[h,h] = de, O[h,h'] = ep this is
[[al, u be], [bp, de + (u - 1) ep]] and de - ep.  The quadratic form of a symmetric O on c_t e_t + c_U 1_U is
c^T Gram [O] c with Gram = diag(1, u) (and (u) on R^U).

Output: the 3 x 3 invariant form Finv on (rho_t, rho_U, q_U) and the 2 x 2 standard form Fstd on (rho_f, q_f) of
Phi~ = D Phi D, as rational functions of (p, q, x) with u = 2 + x, pickled to logs/equ_sym.pkl (numerator and
denominator dictionaries).  Phi is positive definite iff Finv and Fstd are (for u >= 2), and Phi > 0 is equivalent to
the block being >= 0 with kernel exactly p (x) |2@t> (tlib.py)."""
import sys, time, pickle
from flintfield import Field

t0 = time.time()
K = Field(['p', 'q', 'x'])
p, q, x = K.gens
one = K(1)
zero = K(0)
u = x + 2
L = p + u * q
sig = p / L
s = q / L


# ---------------- tiny matrices over K ----------------
def mm(A, B):
    return [[sum((A[i][l] * B[l][j] for l in range(len(B))), zero) for j in range(len(B[0]))] for i in range(len(A))]


def madd(A, B, c=1):
    return [[a + b * c for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mT(A):
    return [list(r) for r in zip(*A)]


def minv(A):
    if len(A) == 1:
        return [[one / A[0][0]]]
    a, b, c, d = A[0][0], A[0][1], A[1][0], A[1][1]
    det = a * d - b * c
    return [[d / det, -b / det], [-c / det, a / det]]


GR = {'k': [[one, zero], [zero, u]], 'U': [[u]]}
GRI = {'k': [[one, zero], [zero, one / u]], 'U': [[one / u]]}


class Op:
    """S_u-invariant operator dom -> cod ('k' = R^[k], 'U' = R^U): matrix M (cod x dom basis coefficients) and the
    standard scalar sd."""

    def __init__(self, dom, cod, M, sd):
        self.dom, self.cod, self.M, self.sd = dom, cod, M, sd

    @staticmethod
    def kk(al, be, bp, de, ep):
        return Op('k', 'k', [[al, u * be], [bp, de + (u - 1) * ep]], de - ep)

    @staticmethod
    def kU(eta, de, ep):          # rows in [k], columns in U; eta = O[t,h]
        return Op('U', 'k', [[u * eta], [de + (u - 1) * ep]], de - ep)

    @staticmethod
    def UU(de, ep):
        return Op('U', 'U', [[de + (u - 1) * ep]], de - ep)

    def __add__(self, o):
        assert (self.dom, self.cod) == (o.dom, o.cod)
        return Op(self.dom, self.cod, madd(self.M, o.M), self.sd + o.sd)

    def __sub__(self, o):
        assert (self.dom, self.cod) == (o.dom, o.cod)
        return Op(self.dom, self.cod, madd(self.M, o.M, -1), self.sd - o.sd)

    def scale(self, c):
        return Op(self.dom, self.cod, [[c * v for v in r] for r in self.M], c * self.sd)

    def __matmul__(self, o):      # plain product
        assert self.dom == o.cod
        return Op(o.dom, self.cod, mm(self.M, o.M), self.sd * o.sd)

    def T(self):
        return Op(self.cod, self.dom, mm(mm(GRI[self.dom], mT(self.M)), GR[self.cod]), self.sd)

    def inv(self):
        assert self.dom == self.cod
        return Op(self.dom, self.cod, minv(self.M), one / self.sd)


# D^2 and D^{-2} (D = diag(w_t, w_h, ...) on R^[k], om I on R^U), in the normalised weights
D2 = {'k': Op.kk(sig, zero, zero, s, zero), 'U': Op.UU(s, zero)}
Dm2 = {'k': Op.kk(one / sig, zero, zero, one / s, zero), 'U': Op.UU(one / s, zero)}


def comp(*ops):
    """product of rescaled operators: O1~ D^{-2} O2~ D^{-2} ..."""
    out = ops[0]
    for o in ops[1:]:
        out = out @ Dm2[out.dom] @ o
    return out


def rinv(o):
    return D2[o.cod] @ o.inv() @ D2[o.dom]


# ---------------- rescaled closed forms (tlib.py), equal U-weights ----------------
# scalars
g_t = sig / (2 * (1 - sig))
gam = g_t + u * s / (2 * (1 - s))
tau = (1 + g_t) / (1 + gam / 2)
gamp = u * s / (2 - sig - 2 * s)
taup = one / (1 + gamp / 2)
# D_R = diag(N): rescaled diag(N_c x_c)
DR = Op.kk(sig, zero, zero, 2 * s, zero)
# Y_RP[c,r] = tau_c/(2(1 - x_c)) (delta_rc - w_r w_c); rescaled tau_c/(2(1-x_c)) (delta x_c - x_r x_c)
ct = (tau - 2) / (2 * (1 - sig))
cu = tau / (2 * (1 - s))
YRP = Op.kk(ct * (sig - sig * sig), ct * (-s * sig), cu * (-sig * s), cu * (s - s * s), cu * (-s * s))
# N^{-1} on the P-states: rescaled diag(x_r/N_r)
Ninv = Op.kk(sig, zero, zero, s / 2, zero)
# Y_PQ[r,c] = taup/(2 - sig - 2x_c) (N_r w_r w_c - (2 - sig) delta_rc); rescaled with x's
cq = taup / (2 - sig - 2 * s)
YPQ = Op.kU(cq * (sig * s), cq * (2 * s * s - (2 - sig) * s), cq * (2 * s * s))
# N1 (k x k): N1[c,c] = sum_{u != c} x_u/(2(1 - x_c - x_u)), N1[c,c'] = w_c w_c'/(2(1 - x_c - x_c'))
N1 = Op.kk(sig * u * s / (2 * (1 - sig - s)), sig * s / (2 * (1 - sig - s)), sig * s / (2 * (1 - sig - s)),
           s * (sig / (2 * (1 - s - sig)) + (u - 1) * s / (2 * (1 - 2 * s))), s * s / (2 * (1 - 2 * s)))
A1 = D2['k'] + N1.scale(K(1) / 2)          # (I + N1/2)~ = D^2 + N1~/2
# N' (U x U): N'[h,h] = sum_{u in U, u != h} x_u/(2 - sig - 2x_h - 2x_u), N'[h,h'] = w w/(2 - sig - 2x_h - 2x_h')
Np = Op.UU(s * (u - 1) * s / (2 - sig - 4 * s), s * s / (2 - sig - 4 * s))
A2 = D2['U'] + Np.scale(K(1) / 2)
# H (k x U): H[h,h] = 1 + sig/(2(1 - sig - x_h)), H[t,h] = w_t w_h/(2(1 - sig - x_h))
H = Op.kU(sig * s / (2 * (1 - sig - s)), s * (1 + sig / (2 * (1 - sig - s))), zero)
YRQ = comp(rinv(A1), H).scale(K(-1))
YQQ = rinv(A2).scale(K(2))

# ---------------- Phi~ = C~ - B~^T N^{-1}~ B~ ----------------
PRR = DR - comp(YRP, Ninv, YRP.T())
PRQ = YRQ - comp(YRP, Ninv, YPQ)
PQQ = YQQ - comp(YPQ.T(), Ninv, YPQ)
PQR = PRQ.T()
# invariant form: Gram * [Phi~] on (rho_t, rho_U, q_U)
Mfull = [PRR.M[0] + PRQ.M[0], PRR.M[1] + PRQ.M[1], PQR.M[0] + PQQ.M[0]]
Gfull = [[one, zero, zero], [zero, u, zero], [zero, zero, u]]
Finv = mm(Gfull, Mfull)
Fstd = [[PRR.sd, PRQ.sd], [PQR.sd, PQQ.sd]]
sym_ok = all(Finv[i][j] == Finv[j][i] for i in range(3) for j in range(3)) and Fstd[0][1] == Fstd[1][0]
print(f'Finv and Fstd built [{time.time() - t0:.1f}s]; symmetric: {sym_ok}', flush=True)


def dd(f):
    return (f.n.to_dict(), f.d.to_dict())


pickle.dump(dict(names=['p', 'q', 'x'], Finv=[[dd(v) for v in r] for r in Finv], Fstd=[[dd(v) for v in r] for r in Fstd]),
            open('logs/equ_sym.pkl', 'wb'))
for i in range(3):
    for j in range(i, 3):
        f = Finv[i][j]
        print(f'Finv[{i}][{j}]: numerator {len(f.n)} terms deg {f.n.total_degree()}, denominator {len(f.d)} terms')
print('PASS' if sym_ok else 'FAIL')
