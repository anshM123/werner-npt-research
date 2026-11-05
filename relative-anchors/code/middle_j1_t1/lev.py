"""Weight levels: the |T| = 1 block of (R_inf)_1 when the weights of the copies in U take l values, symbolically in the
level sizes.  T = {t}; U = U_1 + ... + U_l with |U_i| = u_i and squared weight q_i on U_i; squared weight p on t.
Normalised weights x_t = sigma = p/L, x_h = s_i = q_i/L (h in U_i), L = p + sum u_i q_i.  The closed forms are those of
tlib.py (Lemma T1.3 of middle_j1/NOTES.md, units 2^{k-1}); the unrescaled block is homogeneous of degree 0.

Rescaled frame: D = diag(mu), mu = w_c on R_c, P_c, Q_c; O~ = D O D (rational in p, q_i); (O1 O2)~ = O1~ D^-2 O2~,
(O^-1)~ = D^2 (O~)^-1 D^2, (O^T)~ = (O~)^T.

Symmetry: S_{u_1} x ... x S_{u_l} permutes the copies inside each level.  An invariant operator on R^[k] (or R^U) is
stored by its matrix in the basis (e_t, 1_{U_1}, ..., 1_{U_l}) (or (1_{U_1}, ...)) and by its scalars sd[i] on the
standard parts {f on U_i : sum f = 0}.  Entries of a k x k invariant O: al = O[t,t], be[i] = O[t,h], bp[i] = O[h,t],
de[i] = O[h,h], ep[i][j] = O[h,h'] (h in U_i, h' in U_j, h != h').  Gram matrices diag(1, u_1, ..., u_l) and
diag(u_1, ..., u_l).

build(l, offsets, bscale) returns
  Finv, Fstd : the invariant form of Phi~ on (rho_t, rho_1..rho_l, q_1..q_l) and the standard 2 x 2 forms (rho_f, q_f)
               per level, where Phi = C - B^T diag(1/N) B is the P-first reduction (tlib.py; only for bscale = 1);
  Ginv, Gstd : the same for the gauge-fixed block (the state P_t deleted) on (R, P_U, Q): invariant
               (rho_t, rho_i, pi_i, q_i) and standard 3 x 3 (rho_f, pi_f, q_f) per level, with the B-part of every
               entry multiplied by bscale (A-parts: D_R; diag(N) in Pi; (1/2 - [r = t]) delta in Y_RP; -delta in Y_PQ
               and Y_RQ; 2 I in Y_QQ).
Variables: p, q_1..q_l, x_1..x_l with u_i = offsets[i] + x_i."""
from fractions import Fraction as Fr
import flint
from flintfield import Field


def build(l, offsets, bscale=1, want_phi=True, fixed=None, p_one=False, want_phic=False, check_cong=False):
    """fixed[i] = True: u_i = offsets[i] exactly (no variable x_i).  p_one: dehomogenise p = 1 (every entry is
    homogeneous of degree 0 in (p, q)).  want_phic: also the inverse-free congruent form
    Phi' = diag(A1, A2) Phi diag(A1, A2), A1 = I + N1/2, A2 = I + N'/2 (invariant form 'Cinv', standard 'Cstd')."""
    fixed = fixed or [False] * l
    names = ['p'] + [f'q{i + 1}' for i in range(l)] + [f'x{i + 1}' for i in range(l) if not fixed[i]]
    K = Field(names)
    g = K.gens
    one, zero = K(1), K(0)
    p = one if p_one else g[0]
    qs = g[1:1 + l]
    xg = iter(g[1 + l:])
    us = [K(offsets[i]) if fixed[i] else next(xg) + offsets[i] for i in range(l)]
    L = p + sum((us[i] * qs[i] for i in range(l)), zero)
    sig = p / L
    s = [qq / L for qq in qs]
    def kf(v):
        v = Fr(v)
        return K(flint.fmpq(v.numerator, v.denominator))

    b = kf(bscale)
    LV = range(l)

    # ---------- tiny matrices ----------
    def mm(A, B):
        return [[sum((A[i][m] * B[m][j] for m in range(len(B))), zero) for j in range(len(B[0]))] for i in range(len(A))]

    def madd(A, B, c=1):
        return [[a + bb * c for a, bb in zip(ra, rb)] for ra, rb in zip(A, B)]

    def mT(A):
        return [list(r) for r in zip(*A)]

    def minv(A):
        n = len(A)
        M = [list(A[i]) + [one if i == j else zero for j in range(n)] for i in range(n)]
        for c in range(n):
            r = next(i for i in range(c, n) if M[i][c] != 0)
            M[c], M[r] = M[r], M[c]
            pv = M[c][c]
            M[c] = [v / pv for v in M[c]]
            for i in range(n):
                if i != c and M[i][c] != 0:
                    f = M[i][c]
                    M[i] = [a - f * bb for a, bb in zip(M[i], M[c])]
        return [row[n:] for row in M]

    GR = {'k': [[(one if i == 0 else us[i - 1]) if i == j else zero for j in range(l + 1)] for i in range(l + 1)],
          'U': [[us[i] if i == j else zero for j in range(l)] for i in range(l)]}
    GRI = {'k': [[(one if i == 0 else one / us[i - 1]) if i == j else zero for j in range(l + 1)] for i in range(l + 1)],
           'U': [[one / us[i] if i == j else zero for j in range(l)] for i in range(l)]}

    class Op:
        def __init__(self, dom, cod, M, sd):
            self.dom, self.cod, self.M, self.sd = dom, cod, M, sd

        def __add__(self, o):
            return Op(self.dom, self.cod, madd(self.M, o.M), [a + c for a, c in zip(self.sd, o.sd)])

        def __sub__(self, o):
            return Op(self.dom, self.cod, madd(self.M, o.M, -1), [a - c for a, c in zip(self.sd, o.sd)])

        def scale(self, c):
            return Op(self.dom, self.cod, [[c * v for v in r] for r in self.M], [c * v for v in self.sd])

        def __matmul__(self, o):
            assert self.dom == o.cod
            return Op(o.dom, self.cod, mm(self.M, o.M), [a * c for a, c in zip(self.sd, o.sd)])

        def T(self):
            return Op(self.cod, self.dom, mm(mm(GRI[self.dom], mT(self.M)), GR[self.cod]), list(self.sd))

        def inv(self):
            return Op(self.dom, self.cod, minv(self.M), [one / v for v in self.sd])

        def drop_t_col(self):     # restrict the domain R^[k] -> R^U
            assert self.dom == 'k'
            return Op('U', self.cod, [r[1:] for r in self.M], list(self.sd))

        def drop_t_row(self):     # restrict the codomain R^[k] -> R^U
            assert self.cod == 'k'
            return Op(self.dom, 'U', self.M[1:], list(self.sd))

    def kk(al, be, bp, de, ep):
        M = [[al] + [us[j] * be[j] for j in LV]]
        for i in LV:
            M.append([bp[i]] + [(de[i] + (us[i] - 1) * ep[i][i]) if i == j else us[j] * ep[i][j] for j in LV])
        return Op('k', 'k', M, [de[i] - ep[i][i] for i in LV])

    def kU(eta, de, ep):
        M = [[us[j] * eta[j] for j in LV]]
        for i in LV:
            M.append([(de[i] + (us[i] - 1) * ep[i][i]) if i == j else us[j] * ep[i][j] for j in LV])
        return Op('U', 'k', M, [de[i] - ep[i][i] for i in LV])

    def UU(de, ep):
        M = [[(de[i] + (us[i] - 1) * ep[i][i]) if i == j else us[j] * ep[i][j] for j in LV] for i in LV]
        return Op('U', 'U', M, [de[i] - ep[i][i] for i in LV])

    Z = [zero] * l
    ZZ = [[zero] * l for _ in LV]
    D2 = {'k': kk(sig, Z, Z, list(s), ZZ), 'U': UU(list(s), ZZ)}
    Dm2 = {'k': kk(one / sig, Z, Z, [one / v for v in s], ZZ), 'U': UU([one / v for v in s], ZZ)}

    def comp(*ops):
        out = ops[0]
        for o in ops[1:]:
            out = out @ Dm2[out.dom] @ o
        return out

    def rinv(o):
        return D2[o.cod] @ o.inv() @ D2[o.dom]

    # ---------- rescaled closed forms ----------
    g_t = sig / (2 * (1 - sig))
    gam = g_t + sum((us[i] * s[i] / (2 * (1 - s[i])) for i in LV), zero)
    tau = (1 + g_t) / (1 + gam / 2)
    gamp = sum((us[i] * s[i] / (2 - sig - 2 * s[i]) for i in LV), zero)
    taup = one / (1 + gamp / 2)
    DR = kk(sig, Z, Z, [2 * v for v in s], ZZ)
    ct = (tau - 2) / (2 * (1 - sig))
    cu = [tau / (2 * (1 - s[i])) for i in LV]
    YRP = kk(ct * (sig - sig * sig), [ct * (-s[j] * sig) for j in LV], [cu[i] * (-sig * s[i]) for i in LV],
             [cu[i] * (s[i] - s[i] * s[i]) for i in LV], [[cu[i] * (-s[i] * s[j]) for j in LV] for i in LV])
    Ninv = kk(sig, Z, Z, [v / 2 for v in s], ZZ)
    cq = [taup / (2 - sig - 2 * s[j]) for j in LV]       # column level j
    YPQ = kU([cq[j] * (sig * s[j]) for j in LV], [cq[i] * (2 * s[i] * s[i] - (2 - sig) * s[i]) for i in LV],
             [[cq[j] * (2 * s[i] * s[j]) for j in LV] for i in LV])
    PI = kk(sig - sig * sig / (2 - sig), [-2 * sig * s[j] / (2 - sig) for j in LV],
            [-2 * sig * s[i] / (2 - sig) for i in LV], [2 * s[i] - 4 * s[i] * s[i] / (2 - sig) for i in LV],
            [[-4 * s[i] * s[j] / (2 - sig) for j in LV] for i in LV])
    # N1 (k x k)
    n1tt = sig * sum((us[i] * s[i] / (2 * (1 - sig - s[i])) for i in LV), zero)
    n1th = [sig * s[j] / (2 * (1 - sig - s[j])) for j in LV]
    n1hh = [s[i] * (sig / (2 * (1 - s[i] - sig)) + (us[i] - 1) * s[i] / (2 * (1 - 2 * s[i]))
                    + sum((us[j] * s[j] / (2 * (1 - s[i] - s[j])) for j in LV if j != i), zero)) for i in LV]
    n1hh2 = [[s[i] * s[j] / (2 * (1 - s[i] - s[j])) for j in LV] for i in LV]
    N1 = kk(n1tt, n1th, n1th, n1hh, n1hh2)
    A1 = D2['k'] + N1.scale(kf(Fr(1, 2)))
    nphh = [s[i] * ((us[i] - 1) * s[i] / (2 - sig - 4 * s[i])
                    + sum((us[j] * s[j] / (2 - sig - 2 * s[i] - 2 * s[j]) for j in LV if j != i), zero)) for i in LV]
    nphh2 = [[s[i] * s[j] / (2 - sig - 2 * s[i] - 2 * s[j]) for j in LV] for i in LV]
    Np = UU(nphh, nphh2)
    A2 = D2['U'] + Np.scale(kf(Fr(1, 2)))
    H = kU([sig * s[j] / (2 * (1 - sig - s[j])) for j in LV], [s[i] * (1 + sig / (2 * (1 - sig - s[i]))) for i in LV], ZZ)
    YRQ = comp(rinv(A1), H).scale(K(-1))
    YQQ = rinv(A2).scale(K(2))
    out = dict(K=K, names=names)

    def gram_form(blocks, doms):
        """blocks[a][b]: Op from space doms[b] to doms[a]; returns the invariant form matrix and standard forms."""
        rows = []
        for a, da in enumerate(doms):
            for i in range(len(GR[da])):
                row = []
                for bb, db in enumerate(doms):
                    Mab = mm(GR[da], blocks[a][bb].M)
                    row += Mab[i]
                rows.append(row)
        std = [[[blocks[a][bb].sd[i] for bb in range(len(doms))] for a in range(len(doms))] for i in LV]
        return rows, std

    if want_phic:
        # Phi' = diag(A1, A2) Phi diag(A1, A2) without inverses: A1 Y_RQ = -H and A2 Y_QQ A2 = 2 A2
        CRR = comp(A1, DR, A1) - comp(A1, YRP, Ninv, YRP.T(), A1)
        CRQ = H.scale(K(-1)) @ Dm2['U'] @ A2 - comp(A1, YRP, Ninv, YPQ, A2)
        CQQ = A2.scale(K(2)) - comp(A2, YPQ.T(), Ninv, YPQ, A2)
        out['Cinv'], out['Cstd'] = gram_form([[CRR, CRQ], [CRQ.T(), CQQ]], ['k', 'U'])
    if want_phi:
        PRR = DR - comp(YRP, Ninv, YRP.T())
        PRQ = YRQ - comp(YRP, Ninv, YPQ)
        PQQ = YQQ - comp(YPQ.T(), Ninv, YPQ)
        out['Finv'], out['Fstd'] = gram_form([[PRR, PRQ], [PRQ.T(), PQQ]], ['k', 'U'])
        if want_phic and check_cong:
            def same(X, Y):
                return all(a == c for ra, rc in zip(X.M, Y.M) for a, c in zip(ra, rc)) and all(
                    a == c for a, c in zip(X.sd, Y.sd))
            out['cong_ok'] = (same(comp(A1, PRR, A1), CRR) and same(comp(A1, PRQ, A2), CRQ)
                              and same(comp(A2, PQQ, A2), CQQ))
    # gauge-fixed block with the B-part scaled
    I_k = D2['k']
    I_kU = kU(Z, list(s), ZZ)        # rescaled embedding R^U -> R^[k] (delta), entries x_h on the diagonal
    A_RP = kk(-sig / 2, Z, Z, [v / 2 for v in s], ZZ)
    A_PI = kk(sig, Z, Z, [2 * v for v in s], ZZ)

    def bs(X, A):
        return A + (X - A).scale(b)

    YRPb = bs(YRP, A_RP)
    PIb = bs(PI, A_PI)
    YPQb = bs(YPQ, I_kU.scale(K(-1)))
    YRQb = bs(YRQ, I_kU.scale(K(-1)))
    YQQb = bs(YQQ, D2['U'].scale(K(2)))
    YRPg = YRPb.drop_t_col()                     # k x U
    PIg = PIb.drop_t_col().drop_t_row()          # U x U
    YPQg = YPQb.drop_t_row()                     # U x U
    out['Ginv'], out['Gstd'] = gram_form([[DR, YRPg, YRQb], [YRPg.T(), PIg, YPQg], [YRQb.T(), YPQg.T(), YQQb]],
                                         ['k', 'U', 'U'])
    out['blocks'] = dict(DR=DR, YRP=YRPb, PI=PIb, YPQ=YPQb, YRQ=YRQb, YQQ=YQQb)
    return out


def det(A, zero):
    """determinant by fraction-free-free Gaussian elimination over a field (small matrices)."""
    A = [list(r) for r in A]
    n = len(A)
    d = zero + 1
    for c in range(n):
        r = next((i for i in range(c, n) if A[i][c] != 0), None)
        if r is None:
            return zero
        if r != c:
            A[c], A[r] = A[r], A[c]
            d = -d
        pv = A[c][c]
        d = d * pv
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / pv
                A[i] = [a - f * bb for a, bb in zip(A[i], A[c])]
    return d


def signcheck(f):
    """(sign, description) if f = c prod g_i^e_i with every irreducible g_i one-signed in coefficients and with a
    monomial free of the x-variables (so g_i != 0 on p, q > 0, x >= 0); sign is then the sign of f there.
    (None, reason) otherwise."""
    out = 1
    desc = []
    for poly in (f.n, f.d):
        names = poly.context().names()
        xi = [i for i, nm in enumerate(names) if nm.startswith('x')]
        cs = list(poly.coeffs())
        if (all(v > 0 for v in cs) or all(v < 0 for v in cs)) and any(all(m[i] == 0 for i in xi) for m in poly.monoms()):
            # one-signed as it stands (no factorisation needed)
            out *= 1 if cs[0] > 0 else -1
            desc.append(f'[one-signed {poly.total_degree()}/{len(poly)}]')
            continue
        c, fac = poly.factor()
        sg = 1 if c > 0 else -1
        names = poly.context().names()
        xi = [i for i, nm in enumerate(names) if nm.startswith('x')]
        for gg, e in fac:
            cs = list(gg.coeffs())
            if all(v > 0 for v in cs):
                pass
            elif all(v < 0 for v in cs):
                sg *= (-1) ** e
            else:
                return None, f'factor with mixed signs: deg {gg.total_degree()}, {len(gg)} terms'
            if not any(all(m[i] == 0 for i in xi) for m in gg.monoms()):
                return None, f'factor {gg} vanishes at x = 0'
            desc.append(f'{gg.total_degree()}/{len(gg)}^{e}')
        out *= sg
    return out, ' '.join(desc)
