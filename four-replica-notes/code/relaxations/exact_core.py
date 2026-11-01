"""Exact (rational) two-letter relaxation data at n = 2 for the pilot certificate.

Notation as in PILOT_NOTES.md.  All matrices are python-flint fmpq_mat.  Two-site blocks with row index (i1, i2) and
column index (j1, j2) are vectorised in the order (i1, j1, i2, j2), so that two-site maps are Kronecker products of
per-site maps T[(k, l), (i, j)].
"""
import itertools, numpy as np, flint
from alg4 import S4, PERM, pt, ptrace_last, perm_op, mul, inv
from exact_site import ExactAlg, Q, build_exact, nullspace_q, weights
from exact_maps import exact_A3
from rep4 import realign1

FAMS = ['uuuu', 'uuuv', 'uuvv', 'uvvv', 'vvvv']
fq = flint.fmpq

def locate(ell):
    c = ''.join(sorted(ell))
    for s in S4:
        if all(ell[s[i]] == c[i] for i in range(4)): return c, s
    raise ValueError(ell)

def letter_group(c):
    return [s for s in S4 if all(c[s[i]] == c[i] for i in range(4))]

# ----------------------------------------------------------------------------- small exact helpers
def zeros(r, c): return flint.fmpq_mat(r, c)
def eye(k):
    M = flint.fmpq_mat(k, k)
    for i in range(k): M[i, i] = 1
    return M
def kron(A, B):
    ra, ca, rb, cb = A.nrows(), A.ncols(), B.nrows(), B.ncols()
    M = flint.fmpq_mat(ra * rb, ca * cb)
    for i in range(ra):
        for j in range(ca):
            a = A[i, j]
            if a == 0: continue
            for k in range(rb):
                for l in range(cb):
                    b = B[k, l]
                    if b != 0: M[i * rb + k, j * cb + l] = a * b
    return M
def to_np(M):
    return np.array([[float(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())])
def col(M, j):
    v = flint.fmpq_mat(M.nrows(), 1)
    for i in range(M.nrows()): v[i, 0] = M[i, j]
    return v
def hstack(mats):
    r = mats[0].nrows(); c = sum(m.ncols() for m in mats); M = flint.fmpq_mat(r, c); o = 0
    for m in mats:
        for i in range(r):
            for j in range(m.ncols()): M[i, o + j] = m[i, j]
        o += m.ncols()
    return M
def vstack(mats):
    c = mats[0].ncols(); r = sum(m.nrows() for m in mats); M = flint.fmpq_mat(r, c); o = 0
    for m in mats:
        for i in range(m.nrows()):
            for j in range(c): M[o + i, j] = m[i, j]
        o += m.nrows()
    return M
def rowspace_basis(M):
    """rational basis (as columns) of the row space of M"""
    R, rank = M.rref()
    B = flint.fmpq_mat(M.ncols(), rank)
    for k in range(rank):
        for j in range(M.ncols()): B[j, k] = R[k, j]
    return B
def swap_matrix(f1, f2):
    """P with P e_(i1,i2) = e_(i2,i1): maps the (i1 i2) ordering of C^f1 (x) C^f2 to C^f2 (x) C^f1"""
    P = flint.fmpq_mat(f1 * f2, f1 * f2)
    for i1 in range(f1):
        for i2 in range(f2): P[i2 * f1 + i1, i1 * f2 + i2] = 1
    return P
def vecp(X, f1, f2):
    """vectorise a two-site block (rows (i1,i2), cols (j1,j2)) in the order (i1, j1, i2, j2)"""
    v = flint.fmpq_mat(f1 * f1 * f2 * f2, 1)
    for i1 in range(f1):
        for i2 in range(f2):
            for j1 in range(f1):
                for j2 in range(f2):
                    v[((i1 * f1 + j1) * f2 + i2) * f2 + j2, 0] = X[i1 * f2 + i2, j1 * f2 + j2]
    return v
def unvecp(v, f1, f2, k=0):
    X = flint.fmpq_mat(f1 * f2, f1 * f2)
    for i1 in range(f1):
        for i2 in range(f2):
            for j1 in range(f1):
                for j2 in range(f2):
                    X[i1 * f2 + i2, j1 * f2 + j2] = v[((i1 * f1 + j1) * f2 + i2) * f2 + j2, k]
    return X
def bits(M):
    """max bit size of numerators and denominators of the entries of an fmpq_mat"""
    nb = 0; db = 0
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            q = M[i, j]
            nb = max(nb, int(abs(q.p)).bit_length()); db = max(db, int(q.q).bit_length())
    return nb, db

# ----------------------------------------------------------------------------- per-site exact algebras
def site_algebras():
    A, CB, gensB1, evB = build_exact()
    hwB = {}
    for lam in evB:
        hwB[int(round(lam * 2))] = {3.0: (3, 0, -1), 0.0: (2, 1, -1), -2.0: (2, 0, 0), -4.0: (1, 1, 0)}[round(lam, 6) + 0.0]
    B1 = ExactAlg(gensB1, CB, {3}, hwB, 'B31 exact')
    A3 = exact_A3()
    return A, B1, A3

def exact_coeffs(alg, M):
    """Fourier blocks of an integer matrix M that lies in the algebra (checked exactly)"""
    Mq = Q(M); out = []
    for b in alg.blocks:
        Xh = b['Gi'] * (b['S'].transpose() * (Mq * b['S']))
        assert Mq * b['S'] == b['S'] * Xh, f"{alg.name}: element not in the algebra"
        out.append(Xh)
    return out

class SiteMap:
    """per-site map phi: A -> target algebra, as T[(mu, lam)] = fmpq_mat (g_mu^2 x f_lam^2)"""
    def __init__(self, A, tg, phi):
        C = {s: exact_coeffs(tg, np.rint(phi(PERM[s])).astype(np.int64)) for s in S4}
        self.T = {}
        for li, lb in enumerate(A.blocks):
            f = lb['f']; rinv = {s: lb['rho'][inv(s)] for s in S4}
            for mi, mb in enumerate(tg.blocks):
                g = mb['f']; T = flint.fmpq_mat(g * g, f * f); nz = False
                for i in range(f):
                    for j in range(f):
                        acc = flint.fmpq_mat(g, g)
                        for s in S4:
                            w = rinv[s][j, i]
                            if w != 0: acc = acc + C[s][mi] * w
                        acc = acc * fq(f, 24)
                        for k in range(g):
                            for l in range(g):
                                if acc[k, l] != 0: T[k * g + l, i * f + j] = acc[k, l]; nz = True
                self.T[(mi, li)] = T if nz else None
        # consistency: phi(V_tau) reproduced from rho(tau) for a few words
        for tau in S4[:5]:
            for mi in range(len(tg.blocks)):
                g = tg.blocks[mi]['f']; acc = flint.fmpq_mat(g * g, 1)
                for li, lb in enumerate(A.blocks):
                    if self.T[(mi, li)] is None: continue
                    f = lb['f']; r = lb['rho'][tau]; v = flint.fmpq_mat(f * f, 1)
                    for i in range(f):
                        for j in range(f): v[i * f + j, 0] = r[i, j]
                    acc = acc + self.T[(mi, li)] * v
                Ct = C[tau][mi]
                for k in range(g):
                    for l in range(g): assert acc[k * g + l, 0] == Ct[k, l], "per-site map inconsistent"

# ----------------------------------------------------------------------------- the exact relaxation at n = 2
class ExactRelax2:
    def __init__(self, A, B1, A3, log=print):
        self.A, self.B1, self.A3 = A, B1, A3; nb = len(A.blocks)
        self.f = [b['f'] for b in A.blocks]; self.m = [b['m'] for b in A.blocks]
        self.pats = list(itertools.product(range(nb), repeat=2))
        # ---- variable space
        self.N = 0; self.var = {}         # (c, lam) -> dict(V, Ylist, idx)
        for c in FAMS:
            G = letter_group(c)
            for lam in itertools.combinations_with_replacement(range(nb), 2):
                l1, l2 = lam; f1, f2 = self.f[l1], self.f[l2]; D = f1 * f2
                Pi = zeros(D, D)
                for s in G: Pi = Pi + kron(A.blocks[l1]['rho'][s], A.blocks[l2]['rho'][s])
                Pi = Pi * fq(1, len(G))
                V = rowspace_basis(Pi)             # columns: basis of range(Pi^T)
                t = V.ncols()
                if t == 0: continue
                # symmetric matrices on C^t, with Y = K Y K^T on equal patterns
                syms = []
                for a in range(t):
                    for b_ in range(a, t):
                        E = zeros(t, t); E[a, b_] = 1; E[b_, a] = 1; syms.append(E)
                if l1 == l2:
                    P = swap_matrix(f1, f2)
                    VtV = V.transpose() * V; K = VtV.inv() * (V.transpose() * (P * V))
                    assert P * V == V * K, "swap does not preserve range(Pi^T)"
                    rows = []
                    for E in syms:
                        R = E - K * E * K.transpose()
                        rows.append([R[i, j] for i in range(t) for j in range(t)])
                    Mc = flint.fmpq_mat(len(rows[0]), len(syms), [rows[k][e] for e in range(len(rows[0])) for k in range(len(syms))])
                    Nsp = nullspace_q(Mc)
                    Ys = []
                    for k in range(Nsp.ncols()):
                        Y = zeros(t, t)
                        for e, E in enumerate(syms):
                            if Nsp[e, k] != 0: Y = Y + E * Nsp[e, k]
                        Ys.append(Y)
                else:
                    Ys = syms
                idx = list(range(self.N, self.N + len(Ys))); self.N += len(Ys)
                Gl = kron(A.blocks[l1]['G'], A.blocks[l2]['G']); Gli = kron(A.blocks[l1]['Gi'], A.blocks[l2]['Gi'])
                self.var[(c, lam)] = dict(V=V, Ys=Ys, idx=idx, G=Gl, Gi=Gli, t=t)
        log(f"exact variable space: {self.N} parameters")
        # ---- Xmat[(c, pattern)]: vec'(X_hat_pattern) as a linear function of the parameters (D^2 x N)
        self.Xmat = {}
        for (c, lam), vd in self.var.items():
            l1, l2 = lam; f1, f2 = self.f[l1], self.f[l2]
            cols = []
            for Y in vd['Ys']:
                Xh = vd['Gi'] * (vd['V'] * Y * vd['V'].transpose())
                cols.append(Xh)
            M = zeros(f1 * f1 * f2 * f2, self.N)
            for k, Xh in zip(vd['idx'], cols):
                v = vecp(Xh, f1, f2)
                for e in range(v.nrows()): M[e, k] = v[e, 0]
            self.Xmat[(c, (l1, l2))] = M
            if l1 != l2:
                P = swap_matrix(f1, f2)          # X_hat for the swapped pattern: P X P^T in (i2, i1) ordering
                M2 = zeros(f1 * f1 * f2 * f2, self.N)
                for k, Xh in zip(vd['idx'], cols):
                    v = vecp(P * Xh * P.transpose(), f2, f1)
                    for e in range(v.nrows()): M2[e, k] = v[e, 0]
                self.Xmat[(c, (l2, l1))] = M2
        # ---- functionals
        w01 = PERM[(1, 0, 2, 3)] * 2 - 1 * np.eye(81, dtype=np.int64)      # 2 (F_01 - I/2)
        w12 = PERM[(0, 2, 1, 3)] * 2 - 1 * np.eye(81, dtype=np.int64)
        # q_2 = m_uu + m_vv + 2 m_uv with the factor 2 per site removed (1/4 overall)
        self.c = (self.functional('uuuu', w01, w01) + self.functional('vvvv', w01, w01)
                  + self.functional('uuvv', w12, w12) * 2) * fq(1, 4)
        mrows = []; mrhs = []
        for c in FAMS:
            for tau in S4:
                val = 1; seen = set()
                for i in range(4):
                    if i in seen: continue
                    j = i; letters = set()
                    while j not in seen: seen.add(j); letters.add(c[j]); j = tau[j]
                    if len(letters) > 1: val = 0
                mrows.append(self.functional(c, PERM[tau], PERM[tau])); mrhs.append(val)
        self.mom = (vstack(mrows), flint.fmpq_mat(len(mrhs), 1, mrhs))
        log("functionals done")
        # ---- identities (target A3): ORTH and MARG
        idrows = []
        for c in FAMS:
            for i, j in itertools.combinations(range(4), 2):
                if c[i] == c[j]: continue
                sw = list(range(4)); sw[i], sw[j] = j, i; Fij = PERM[tuple(sw)]
                order = [r for r in range(4) if r != j] + [j]; mv = PERM[tuple(order.index(r) for r in range(4))]
                sm = SiteMap(A, A3, lambda X, F=Fij, M=mv: ptrace_last(M @ (X @ F) @ M.T))
                idrows.append(self.linear_image(A3, [(fq(1), c, sm)]))
        for xyz in ['uuu', 'uuv', 'uvv', 'vvv']:
            terms = []
            for last, sg in (('u', 1), ('v', -1)):
                c, s = locate(xyz + last); P = PERM[s]
                terms.append((fq(sg), c, SiteMap(A, A3, lambda X, P=P: ptrace_last(P @ X @ P.T))))
            idrows.append(self.linear_image(A3, terms))
        Eid = vstack(idrows)
        R, rank = Eid.rref()
        self.ident = flint.fmpq_mat(rank, self.N, [R[i, j] for i in range(rank) for j in range(self.N)])
        log(f"identities: {Eid.nrows()} rows, rank {rank}")
        # ---- PSD blocks: X blocks (Y_lam >= 0) and covariance blocks
        self.psd = []        # list of (name, size, list over parameters of symmetric fmpq_mat)
        for (c, lam), vd in self.var.items():
            t = vd['t']; mats = [zeros(t, t) for _ in range(self.N)]
            for k, Y in zip(vd['idx'], vd['Ys']): mats[k] = Y
            self.psd.append((f"X{c}{lam}", t, mats))
        for xy in ['uu', 'uv', 'vv']:
            maps = {}
            for r, lr in enumerate('uv'):
                for s_, ls in enumerate('uv'):
                    c, s = locate(xy + lr + ls); P = PERM[s]
                    maps[(r, s_)] = (c, SiteMap(A, B1, lambda X, P=P: realign1(P @ X @ P.T)))
            for mu in itertools.combinations_with_replacement(range(len(B1.blocks)), 2):
                g1, g2 = B1.blocks[mu[0]]['f'], B1.blocks[mu[1]]['f']; D = g1 * g2
                Gt = kron(B1.blocks[mu[0]]['G'], B1.blocks[mu[1]]['G'])
                ent = {}
                for rs, (c, sm) in maps.items():
                    ent[rs] = self.target_image(B1, mu, [(fq(1), c, sm)])      # vec' (D^2 x N)
                mats = []
                nz = False
                for k in range(self.N):
                    Bk = zeros(2 * D, 2 * D)
                    for (r, s_), Mv in ent.items():
                        Yh = unvecp(Mv, g1, g2, k)
                        Bl = Gt * Yh
                        for a in range(D):
                            for b_ in range(D): Bk[r * D + a, s_ * D + b_] = Bl[a, b_]
                    Bk = (Bk + Bk.transpose()) * fq(1, 2)
                    if not nz and any(Bk[a, b_] != 0 for a in range(2 * D) for b_ in range(2 * D)): nz = True
                    mats.append(Bk)
                if nz: self.psd.append((f"COV{xy}{mu}", 2 * D, mats))
        log(f"PSD blocks: {len(self.psd)} (max size {max(p[1] for p in self.psd)})")

    def functional(self, c, w1, w2):
        """row vector (1 x N): x -> Tr[X_c (w1 (x) w2)]"""
        h1 = exact_coeffs(self.A, w1); h2 = exact_coeffs(self.A, w2)
        out = zeros(1, self.N)
        for (l1, l2) in self.pats:
            key = (c, (l1, l2))
            if key not in self.Xmat: continue
            f1, f2 = self.f[l1], self.f[l2]; W = kron(h1[l1], h2[l2]); mlam = self.m[l1] * self.m[l2]
            # Tr(X W) = sum_{I,J} X[I,J] W[J,I]
            Wt = W.transpose(); wv = vecp(Wt, f1, f2)
            out = out + (wv.transpose() * self.Xmat[key]) * mlam
        return out

    def target_image(self, tg, mu, terms):
        """vec' of the target block mu of sum coef * phi(X_c), as a (D_mu^2 x N) matrix"""
        g1, g2 = tg.blocks[mu[0]]['f'], tg.blocks[mu[1]]['f']
        out = zeros(g1 * g1 * g2 * g2, self.N)
        for coef, c, sm in terms:
            for (l1, l2) in self.pats:
                key = (c, (l1, l2))
                if key not in self.Xmat: continue
                T1 = sm.T[(mu[0], l1)]; T2 = sm.T[(mu[1], l2)]
                if T1 is None or T2 is None: continue
                out = out + kron(T1, T2) * self.Xmat[key] * coef
        return out

    def linear_image(self, tg, terms):
        rows = [self.target_image(tg, mu, terms) for mu in itertools.product(range(len(tg.blocks)), repeat=2)]
        return vstack(rows)

    # ------------------------------------------------------------------ exact true points
    def site_hat(self, vecs):
        """exact per-site Fourier blocks X_hat_lam of twirl((x)_r |l_r><l_r|) for real rational vectors l_r in Q^3"""
        # word moments m_sigma = Tr[(x)|l_r><l_r| V_sigma] = sum_x prod_r l_r[x_r] l_r[pi_sigma(x)_r]
        out = []
        idx = list(itertools.product(range(3), repeat=4))
        pos = {x: k for k, x in enumerate(idx)}
        msig = {}
        for s in S4:
            tot = fq(0)
            for x in idx:
                y = [0] * 4
                for r in range(4): y[s[r]] = x[r]
                term = fq(1)
                for r in range(4):
                    term *= vecs[r][x[r]] * vecs[r][y[r]]
                    if term == 0: break
                tot += term
            msig[s] = tot
        norm = fq(1)
        for r in range(4): norm *= sum((fq(z) * fq(z) for z in vecs[r]), fq(0))
        for s in S4: msig[s] = msig[s] / norm
        for b in self.A.blocks:
            f = b['f']; X = zeros(f, f)
            for s in S4:
                X = X + b['rho'][inv(s)] * msig[s]
            out.append(X * fq(f, 24 * b['m']))
        return out

    def x_from_hats(self, hats):
        """parameters of a family-wise collection of two-site Fourier blocks: hats[c][(l1,l2)] = X_hat (exact);
        returns (x, ok) with ok = exact reconstruction (swap support and site symmetry hold exactly)"""
        x = zeros(self.N, 1); ok = True
        for (c, lam), vd in self.var.items():
            Xh = hats[c][lam]; B = vd['G'] * Xh
            V = vd['V']; VtV = V.transpose() * V; Vp = VtV.inv() * V.transpose()
            Y = Vp * B * Vp.transpose()
            if not (V * Y * V.transpose() == B): ok = False
            # coordinates of Y in the basis Ys (exact least squares on the flattened basis)
            t = vd['t']; Bm = flint.fmpq_mat(t * t, len(vd['Ys']), [vd['Ys'][k][i, j] for i in range(t) for j in range(t) for k in range(len(vd['Ys']))])
            yv = flint.fmpq_mat(t * t, 1, [Y[i, j] for i in range(t) for j in range(t)])
            G = Bm.transpose() * Bm; coef = G.inv() * (Bm.transpose() * yv)
            if not (Bm * coef == yv): ok = False
            for k, i in enumerate(vd['idx']): x[i, 0] = coef[k, 0]
        # patterns without variables must carry no weight
        for c in FAMS:
            for lam in itertools.combinations_with_replacement(range(len(self.A.blocks)), 2):
                if (c, lam) not in self.var and hats[c][lam] != zeros(hats[c][lam].nrows(), hats[c][lam].ncols()): ok = False
        return x, ok

    def evaluate(self, x):
        """exact checks at a point: moment residual, identity residual, objective, smallest block eigenvalue (float)"""
        M, e = self.mom
        mres = M * x - e; ires = self.ident * x
        obj = (self.c * x)[0, 0]
        mins = {}
        for name, size, mats in self.psd:
            B = zeros(size, size)
            for k in range(self.N):
                if x[k, 0] != 0: B = B + mats[k] * x[k, 0]
            ev = np.linalg.eigvalsh(to_np(B)).min()
            mins[name[:3]] = min(mins.get(name[:3], np.inf), ev)
        return mres, ires, obj, mins


def two_site_hats(rel, Ls):
    """exact site-symmetrised two-site Fourier blocks of twirl((x)_r |l_r><l_r|) for real integer vectors l_r in Z^9
    (given as 3x3 integer matrices L_r[x1, x2]); returns dict pattern (all 16) -> X_hat"""
    A = rel.A; nb = len(A.blocks)
    norm = fq(1)
    for L in Ls: norm *= fq(int(np.sum(L.astype(object) * L.astype(object))))
    letters = 'abcdefgh'
    # word moments m[s1, s2] = sum_x prod_r L_r[x^r_1, x^r_2] L_r[y^r_1, y^r_2],  y^{s(r)} = x^r per site
    m = {}
    for s1 in S4:
        for s2 in S4:
            # labels: x^r_site -> letter; y^r_site = x^{s^-1(r)}_site
            subs = []
            for r in range(4):
                subs.append(letters[r] + letters[4 + r])
                subs.append(letters[inv(s1)[r]] + letters[4 + inv(s2)[r]])
            ops = []
            for r in range(4): ops += [Ls[r].astype(object), Ls[r].astype(object)]      # exact Python integers
            m[(s1, s2)] = int(np.einsum(','.join(subs) + '->', *ops))
    R = {}
    for li, b in enumerate(A.blocks):
        f = b['f']
        for s in S4: R[(li, s)] = b['rho'][inv(s)] * fq(f, 24 * b['m'])
    hats = {}
    for l1 in range(nb):
        for l2 in range(nb):
            f1, f2 = A.blocks[l1]['f'], A.blocks[l2]['f']; X = zeros(f1 * f2, f1 * f2)
            for (s1, s2), val in m.items():
                if val == 0: continue
                X = X + kron(R[(l1, s1)], R[(l2, s2)]) * fq(val)
            hats[(l1, l2)] = X * (fq(1) / norm)
    sym = {}
    for (l1, l2), X in hats.items():
        f1, f2 = A.blocks[l1]['f'], A.blocks[l2]['f']; P = swap_matrix(f2, f1)
        sym[(l1, l2)] = (X + P * hats[(l2, l1)] * P.transpose()) * fq(1, 2)
    return sym

def exact_q2(P):
    """exact q_2(P) for a rational 9x9 matrix (two sites of dimension 3)"""
    from fractions import Fraction
    T = P.reshape(3, 3, 3, 3)
    tot = Fraction(0)
    fro = lambda M: sum((Fraction(z) * Fraction(z) for z in np.asarray(M, dtype=object).ravel()), Fraction(0))
    tot += fro(P)
    tot -= Fraction(1, 2) * fro(np.einsum('abcb->ac', T)) + Fraction(1, 2) * fro(np.einsum('abad->bd', T))
    tot += Fraction(1, 4) * Fraction(np.einsum('abab->', T)) ** 2
    return tot
