"""Exact two-letter relaxation data for n sites (n = 2, 3), for the certificate pipeline.

Same conventions as exact_core.py (PILOT_NOTES.md), generalised to n sites:
  variables   per family c and sorted site pattern lam: B_lam = V Y V^T (V exact basis of range(Pi^T), Y symmetric and
              invariant under the site stabiliser of lam);  X_hat_lam = G_lam^{-1} B_lam;
  tensors     the coefficients of X_hat are stored as object arrays (gmpy2.mpq) with axes (K, i0, j0, i1, j1, ...),
              so that per-site maps act on the pairs (i_s, j_s) as matrices (mode products);
  maps        per-site maps from SiteMap (exact_core) applied site by site.
"""
import itertools, time, numpy as np, gmpy2, flint
from gmpy2 import mpq
from alg4 import S4, PERM, ptrace_last, mul, inv
from exact_site import nullspace_q
from exact_core import (site_algebras, exact_coeffs, SiteMap, FAMS, locate, letter_group, rowspace_basis, fq)
from rep4 import realign1

def f2m(q): return mpq(int(q.p), int(q.q))
def fm2np(M): return np.array([[f2m(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())], dtype=object).reshape(M.nrows(), M.ncols())
def np2fm(A):
    A = np.asarray(A, dtype=object)
    return flint.fmpq_mat(A.shape[0], A.shape[1], [flint.fmpq(int(gmpy2.mpq(x).numerator), int(gmpy2.mpq(x).denominator)) for x in A.ravel()])
def okron(mats):
    M = np.array([[mpq(1)]], dtype=object)
    for X in mats: M = np.kron(M, X)
    return M
def zeros_o(shape):
    Z = np.empty(shape, dtype=object); Z.fill(mpq(0)); return Z
def perm_matrix_o(shape, p):
    D = int(np.prod(shape)); E = np.zeros((D,) + tuple(shape), dtype=np.int64)
    for k in range(D): E.reshape(D, D)[k, k] = 1
    E = E.reshape((D,) + tuple(shape)).transpose([0] + [1 + a for a in p]).reshape(D, D)
    return E.T.astype(object) * mpq(1)

class ExactRelaxN:
    def __init__(self, n, A, B1, A3, log=print):
        t0 = time.time()
        self.n = n; self.A, self.B1, self.A3 = A, B1, A3; nb = len(A.blocks)
        self.f = [b['f'] for b in A.blocks]; self.m = [b['m'] for b in A.blocks]
        self.rho = {(l, s): fm2np(A.blocks[l]['rho'][s]) for l in range(nb) for s in S4}
        self.Gs = [fm2np(A.blocks[l]['G']) for l in range(nb)]; self.Gis = [fm2np(A.blocks[l]['Gi']) for l in range(nb)]
        self.N = 0; self.var = {}; self.src = {c: [] for c in FAMS}
        for c in FAMS:
            Gc = letter_group(c)
            for lam in itertools.combinations_with_replacement(range(nb), n):
                shape = tuple(self.f[l] for l in lam); D = int(np.prod(shape))
                Pi = sum((okron([self.rho[(l, s)] for l in lam]) for s in Gc), zeros_o((D, D))) * mpq(1, len(Gc))
                V = fm2np(rowspace_basis(np2fm(Pi))); t = V.shape[1]
                if t == 0: continue
                stab = [p for p in itertools.permutations(range(n)) if all(lam[p[q]] == lam[q] for q in range(n))]
                Vq = np2fm(V); VtVi = (Vq.transpose() * Vq).inv()
                Ks = []
                for p in stab:
                    if list(p) == list(range(n)): continue
                    Pq = np2fm(perm_matrix_o(shape, p))
                    K = VtVi * (Vq.transpose() * (Pq * Vq))
                    assert Pq * Vq == Vq * K, "site permutation does not preserve range(Pi^T)"
                    Ks.append(K)
                syms = []
                for a in range(t):
                    for b_ in range(a, t):
                        E = flint.fmpq_mat(t, t); E[a, b_] = 1; E[b_, a] = 1; syms.append(E)
                if Ks:
                    rows = []
                    for E in syms:
                        r = []
                        for K in Ks:
                            Rm = E - K * E * K.transpose(); r += [Rm[i, j] for i in range(t) for j in range(t)]
                        rows.append(r)
                    nr = len(rows[0])
                    Mc = flint.fmpq_mat(nr, len(syms), [rows[k][e] for e in range(nr) for k in range(len(syms))])
                    Nsp = nullspace_q(Mc)
                    Ys = []
                    for k in range(Nsp.ncols()):
                        Y = flint.fmpq_mat(t, t)
                        for e, E in enumerate(syms):
                            if Nsp[e, k] != 0: Y = Y + E * Nsp[e, k]
                        Ys.append(Y)
                else:
                    Ys = syms
                if not Ys: continue
                idx = list(range(self.N, self.N + len(Ys))); self.N += len(Ys)
                Gil = okron([self.Gis[l] for l in lam]); Gl = okron([self.Gs[l] for l in lam])
                Xh = np.array([Gil.dot(V.dot(fm2np(Y)).dot(V.T)) for Y in Ys], dtype=object)      # (K, D, D)
                self.var[(c, lam)] = dict(V=V, Ys=Ys, idx=idx, G=Gl, Gi=Gil, t=t, Xh=Xh, shape=shape)
                # tensors for every pattern in the site orbit: axes (K, i0, j0, i1, j1, ...)
                Xt = Xh.reshape((len(Ys),) + shape + shape)
                pair_axes = [0] + [a for q in range(n) for a in (1 + q, 1 + n + q)]
                seen = set()
                for perm in itertools.permutations(range(n)):
                    lp = tuple(lam[perm[q]] for q in range(n))
                    if lp in seen: continue
                    seen.add(lp)
                    Tp = Xt.transpose([0] + [1 + perm[q] for q in range(n)] + [1 + n + perm[q] for q in range(n)])
                    self.src[c].append(dict(lam=lp, T=np.ascontiguousarray(Tp.transpose(pair_axes)), idx=idx))
        log(f"exact variable space n={n}: {self.N} parameters; {time.time() - t0:.1f}s")

    # ------------------------------------------------------------------ helpers
    def apply_sites(self, T, mats):
        """T: (K, a0, b0, a1, b1, ...); mats[s]: object array (g, g, f, f) acting on (a_s, b_s) -> (k_s, l_s)"""
        n = self.n; Y = T
        for s in range(n):
            ax = 1 + 2 * s
            Y = np.tensordot(Y, mats[s], axes=([ax, ax + 1], [2, 3]))         # new axes appended at the end
            Y = np.moveaxis(Y, [-2, -1], [ax, ax + 1])
        return Y

    def functional(self, c, ws):
        """row (N,) of x -> Tr[X_c ((x)_s ws[s])], ws: per-site integer matrices in A"""
        h = [exact_coeffs(self.A, w) for w in ws]
        hs = [[fm2np(hl) for hl in hh] for hh in h]
        out = zeros_o((self.N,))
        for sd in self.src[c]:
            lam = sd['lam']; mlam = int(np.prod([self.m[l] for l in lam]))
            mats = []
            for s, l in enumerate(lam):
                W = hs[s][l]; f = W.shape[0]
                # Tr(X W) = sum X[i,j] W[j,i]: per site a (1,1,f,f) map with entries W^T
                mats.append(W.T.reshape(1, 1, f, f))
            Y = self.apply_sites(sd['T'], mats)
            vals = Y.reshape(len(sd['idx']), -1).sum(axis=1)
            for k, i in enumerate(sd['idx']): out[i] += vals[k] * mlam
        return out

    def sitemap_obj(self, sm):
        return {key: (None if T is None else fm2np(T).reshape(int(round(T.nrows() ** 0.5)), int(round(T.nrows() ** 0.5)),
                                                              int(round(T.ncols() ** 0.5)), int(round(T.ncols() ** 0.5))))
                for key, T in sm.T.items()}

    def target_image(self, tg, mu, terms):
        """dict param index -> object array (D_mu, D_mu) of the target Fourier block mu of sum coef * phi(X_c)"""
        n = self.n; g = [tg.blocks[q]['f'] for q in mu]; Dm = int(np.prod(g))
        acc = {}
        for coef, c, smo in terms:
            for sd in self.src[c]:
                lam = sd['lam']
                mats = [smo[(mu[s], lam[s])] for s in range(n)]
                if any(M_ is None for M_ in mats): continue
                Y = self.apply_sites(sd['T'], mats)            # (K, k0, l0, k1, l1, ...)
                K = Y.shape[0]
                Y = Y.reshape((K,) + tuple(x for q in range(n) for x in (g[q], g[q])))
                Y = Y.transpose([0] + [1 + 2 * q for q in range(n)] + [2 + 2 * q for q in range(n)]).reshape(K, Dm, Dm)
                for k, i in enumerate(sd['idx']):
                    acc[i] = acc[i] + Y[k] * coef if i in acc else Y[k] * coef
        return acc

    def build(self, log=print):
        t0 = time.time(); n = self.n; A, B1, A3 = self.A, self.B1, self.A3
        w01 = PERM[(1, 0, 2, 3)] * 2 - np.eye(81, dtype=np.int64); w12 = PERM[(0, 2, 1, 3)] * 2 - np.eye(81, dtype=np.int64)
        sc = mpq(1, 2 ** n)
        self.c = (self.functional('uuuu', [w01] * n) + self.functional('vvvv', [w01] * n) + self.functional('uuvv', [w12] * n) * 2) * sc
        mrows = []; mrhs = []
        for c in FAMS:
            for tau in S4:
                val = 1; seen = set()
                for i in range(4):
                    if i in seen: continue
                    j = i; letters = set()
                    while j not in seen: seen.add(j); letters.add(c[j]); j = tau[j]
                    if len(letters) > 1: val = 0
                mrows.append(self.functional(c, [PERM[tau]] * n)); mrhs.append(val)
        self.mom = (np.array(mrows, dtype=object), np.array([mpq(v) for v in mrhs], dtype=object))
        log(f"   functionals; {time.time() - t0:.1f}s")
        # identities (sorted target patterns of the S_3 algebra)
        idrows = []
        maps = []
        for c in FAMS:
            for i, j in itertools.combinations(range(4), 2):
                if c[i] == c[j]: continue
                sw = list(range(4)); sw[i], sw[j] = j, i; Fij = PERM[tuple(sw)]
                order = [r for r in range(4) if r != j] + [j]; mv = PERM[tuple(order.index(r) for r in range(4))]
                maps.append([(mpq(1), c, self.sitemap_obj(SiteMap(A, A3, lambda X, F=Fij, M=mv: ptrace_last(M @ (X @ F) @ M.T))))])
        for xyz in ['uuu', 'uuv', 'uvv', 'vvv']:
            terms = []
            for last, sg in (('u', 1), ('v', -1)):
                c, s = locate(xyz + last); P = PERM[s]
                terms.append((mpq(sg), c, self.sitemap_obj(SiteMap(A, A3, lambda X, P=P: ptrace_last(P @ X @ P.T)))))
            maps.append(terms)
        for terms in maps:
            for mu in itertools.combinations_with_replacement(range(len(A3.blocks)), n):
                acc = self.target_image(A3, mu, terms)
                if not acc: continue
                Dm = next(iter(acc.values())).shape[0]
                R = zeros_o((Dm * Dm, self.N))
                for i, Y in acc.items(): R[:, i] = Y.ravel()
                idrows.append(R)
        Eid = np.vstack(idrows)
        Fq = np2fm(Eid); Rr, rank = Fq.rref()
        self.ident = np.array([[f2m(Rr[i, j]) for j in range(self.N)] for i in range(rank)], dtype=object).reshape(rank, self.N)
        log(f"   identities: {Eid.shape[0]} rows, rank {rank}; {time.time() - t0:.1f}s")
        # PSD blocks: X blocks and covariance blocks (sorted target patterns, raw)
        self.psd = []      # (name, size, dict param -> object array (size, size))
        for (c, lam), vd in self.var.items():
            d = {}
            for k, i in enumerate(vd['idx']): d[i] = fm2np(vd['Ys'][k])
            self.psd.append((f"X{c}{lam}", vd['t'], d))
        for xy in ['uu', 'uv', 'vv']:
            ent = {}
            for r, lr in enumerate('uv'):
                for s_, ls in enumerate('uv'):
                    c, s = locate(xy + lr + ls); P = PERM[s]
                    ent[(r, s_)] = (c, self.sitemap_obj(SiteMap(A, B1, lambda X, P=P: realign1(P @ X @ P.T))))
            for mu in itertools.combinations_with_replacement(range(len(B1.blocks)), n):
                g = [B1.blocks[q]['f'] for q in mu]; Dm = int(np.prod(g))
                Gt = okron([fm2np(B1.blocks[q]['G']) for q in mu])
                imgs = {rs: self.target_image(B1, mu, [(mpq(1), c, smo)]) for rs, (c, smo) in ent.items()}
                params = set().union(*[set(v.keys()) for v in imgs.values()])
                if not params: continue
                d = {}
                for i in params:
                    Bk = zeros_o((2 * Dm, 2 * Dm))
                    for (r, s_), acc in imgs.items():
                        if i in acc: Bk[r * Dm:(r + 1) * Dm, s_ * Dm:(s_ + 1) * Dm] = Gt.dot(acc[i])
                    Bk = (Bk + Bk.T) * mpq(1, 2)
                    if any(x != 0 for x in Bk.ravel()): d[i] = Bk
                if d: self.psd.append((f"COV{xy}{mu}", 2 * Dm, d))
            log(f"   covariance blocks for {xy}; {time.time() - t0:.1f}s")
        log(f"   PSD blocks: {len(self.psd)} (max size {max(p[1] for p in self.psd)}); {time.time() - t0:.1f}s")

    # ------------------------------------------------------------------ points
    def x_from_hats(self, hats):
        """hats[c][lam] (sorted lam) exact object arrays X_hat (D x D); returns (x, ok)"""
        x = zeros_o((self.N,)); ok = True
        for (c, lam), vd in self.var.items():
            B = vd['G'].dot(hats[c][lam]); V = vd['V']
            Vq = np2fm(V); Vp = (Vq.transpose() * Vq).inv() * Vq.transpose()
            Bq = np2fm(B); Yq = Vp * Bq * Vp.transpose()
            if not (Vq * Yq * Vq.transpose() == Bq): ok = False
            t = vd['t']; Bm = flint.fmpq_mat(t * t, len(vd['Ys']), [vd['Ys'][k][i, j] for i in range(t) for j in range(t) for k in range(len(vd['Ys']))])
            yv = flint.fmpq_mat(t * t, 1, [Yq[i, j] for i in range(t) for j in range(t)])
            coef = (Bm.transpose() * Bm).inv() * (Bm.transpose() * yv)
            if not (Bm * coef == yv): ok = False
            for k, i in enumerate(vd['idx']): x[i] = f2m(coef[k, 0])
        for c in FAMS:
            for lam in itertools.combinations_with_replacement(range(len(self.A.blocks)), self.n):
                if (c, lam) not in self.var and any(v != 0 for v in hats[c][lam].ravel()): ok = False
        return x, ok

    def block_at(self, blk, x):
        name, size, d = blk; B = zeros_o((size, size))
        for i, M in d.items():
            if x[i] != 0: B = B + M * x[i]
        return B

    def evaluate(self, x):
        Mm, em = self.mom
        mres = Mm.dot(x) - em; ires = self.ident.dot(x); obj = self.c.dot(x)
        mins = {}
        for blk in self.psd:
            B = self.block_at(blk, x)
            ev = np.linalg.eigvalsh(np.array(B, dtype=float)).min()
            mins[blk[0][:3]] = min(mins.get(blk[0][:3], np.inf), ev)
        return mres, ires, obj, mins

def fourier_from_moments(rel, m):
    """exact X_hat for all patterns from the word-moment tensor m (object array (24,)*n of mpq):
    X_hat_lam[(i),(j)] = sum_sigma m_sigma prod_s (f/(24 m)) rho_lam_s(sigma_s^{-1})[i_s, j_s]"""
    n = rel.n; nb = len(rel.A.blocks)
    F = []          # per site: list over lam of object array (f, f, 24)
    for l in range(nb):
        f = rel.f[l]; Fl = zeros_o((f, f, 24))
        for k, s in enumerate(S4):
            Fl[:, :, k] = rel.rho[(l, inv(s))] * mpq(f, 24 * rel.m[l])
        F.append(Fl)
    hats = {}
    for lam in itertools.product(range(nb), repeat=n):
        Y = m
        for s in range(n):
            # contract the s-th word axis (currently at position 2*s) with F[lam_s] -> axes (i_s, j_s)
            Y = np.tensordot(Y, F[lam[s]], axes=([2 * s], [2]))       # appended (i, j) at the end
            Y = np.moveaxis(Y, [-2, -1], [2 * s, 2 * s + 1])
        shape = tuple(rel.f[l] for l in lam); D = int(np.prod(shape))
        Y = Y.reshape(tuple(x for q in range(n) for x in (shape[q], shape[q])))
        Y = Y.transpose([2 * q for q in range(n)] + [2 * q + 1 for q in range(n)]).reshape(D, D)
        hats[lam] = Y
    return hats

def symmetrise_sites(rel, hats):
    """site-symmetrised X_hat for sorted patterns from X_hat for all patterns"""
    n = rel.n; out = {}
    perms = list(itertools.permutations(range(n)))
    for lam in itertools.combinations_with_replacement(range(len(rel.A.blocks)), n):
        shape = tuple(rel.f[l] for l in lam); D = int(np.prod(shape)); acc = zeros_o((D, D))
        for p in perms:
            # pattern lam' with lam'[q] = lam[p[q]]: its X_hat transported back to lam ordering
            lp = tuple(lam[p[q]] for q in range(n)); shp = tuple(rel.f[l] for l in lp)
            Y = hats[lp].reshape(shp + shp)
            pinv = [0] * n
            for q in range(n): pinv[p[q]] = q
            Y = Y.transpose([pinv[q] for q in range(n)] + [n + pinv[q] for q in range(n)]).reshape(D, D)
            acc = acc + Y
        out[lam] = acc * mpq(1, len(perms))
    return out


# ----------------------------------------------------------------------------- exact word moments via primes + CRT
PRIMES = [268435399, 268435367, 268435361, 268435331, 268435313, 268435291, 268435273, 268435243, 268435183,
          268435171, 268435157, 268435147, 268435133, 268435129, 268435121, 268435109]

def _contract_mod(ops, subs, out, p):
    """contract a network of int64 tensors modulo p by greedy pairwise einsum, reducing after each step"""
    ops = [o % p for o in ops]; subs = list(subs)
    while len(ops) > 1:
        best = None
        for i in range(len(ops)):
            for j in range(i + 1, len(ops)):
                sh = set(subs[i]) & set(subs[j])
                if sh:
                    others = set(''.join(subs[k] for k in range(len(ops)) if k not in (i, j))) | set(out)
                    keep = ''.join(ch for ch in subs[i] + subs[j] if (ch not in sh or ch in others))
                    keep = ''.join(dict.fromkeys(keep))
                    size = 1
                    for ch in keep: size *= 3
                    if best is None or size < best[0]: best = (size, i, j, keep)
        if best is None:
            i, j = 0, 1; keep = ''.join(dict.fromkeys(subs[0] + subs[1]))
        else:
            _, i, j, keep = best
        res = np.einsum(subs[i] + ',' + subs[j] + '->' + keep, ops[i], ops[j]) % p
        ops = [ops[k] for k in range(len(ops)) if k not in (i, j)] + [res]
        subs = [subs[k] for k in range(len(subs)) if k not in (i, j)] + [keep]
    res = ops[0]
    if subs[0] != out: res = np.einsum(subs[0] + '->' + out, res) % p
    return res

def moments_modp(Ls, n, p):
    """word moments Tr[(x)_r |l_r><l_r| (sigma_1 (x) ... (x) sigma_n)] mod p for integer vectors l_r in Z^(3^n)
    (unnormalised); returns int64 array (24,)*n"""
    letters = [chr(c) for c in list(range(65, 91)) + list(range(97, 123))]
    a = letters[:4]; b = letters[4:8]
    X = [[letters[8 + r * n + s] for s in range(n)] for r in range(4)]
    T = [L.reshape((3,) * n).astype(np.int64) % p for L in Ls]
    PTf = np.array([PERM[s].T.ravel() for s in S4], dtype=np.int64)
    out = np.zeros((24,) * n, dtype=np.int64)
    for head in itertools.product(range(24), repeat=n - 1):
        sig = [S4[k] for k in head]
        subs = []; ops = []
        for r in range(4):
            subs.append(''.join(X[r][s] for s in range(n - 1)) + a[r]); ops.append(T[r])
            subs.append(''.join(X[inv(sig[s])[r]][s] for s in range(n - 1)) + b[r]); ops.append(T[r])
        M = _contract_mod(ops, subs, ''.join(a) + ''.join(b), p).reshape(-1)
        out[head] = (PTf @ M) % p
    return out

def moments_exact(Ls, n):
    """exact normalised word moments (object array of mpq) via several primes and CRT"""
    import math
    bound = 1
    for L in Ls: bound *= int(np.sum(np.abs(L.astype(object)))) ** 2
    prod = 1; res = []; ps = []
    for p in PRIMES:
        res.append(moments_modp(Ls, n, p)); ps.append(p); prod *= p
        if prod > 4 * bound: break
    assert prod > 4 * bound, "not enough primes"
    # CRT to symmetric residues
    shape = res[0].shape; flat = [r.ravel() for r in res]
    vals = np.empty(flat[0].size, dtype=object)
    for idx in range(flat[0].size):
        x = 0; M = 1
        for r, p in zip(flat, ps):
            # combine x mod M with r mod p
            t = ((int(r[idx]) - x) * pow(M, -1, p)) % p
            x = x + M * t; M *= p
        if x > M // 2: x -= M
        vals[idx] = x
    norm = 1
    for L in Ls: norm *= int(np.sum(L.astype(object) * L.astype(object)))
    return np.array([mpq(v, norm) for v in vals], dtype=object).reshape(shape), len(ps)

def exact_qn(P, n):
    """exact q_n(P) for a rational matrix given as an object array of mpq (3^n x 3^n)"""
    T = P.reshape((3,) * (2 * n)); tot = mpq(0)
    for k in range(n + 1):
        for J in itertools.combinations(range(n), k):
            rows = [chr(97 + s) for s in range(n)]; cols = [chr(97 + n + s) for s in range(n)]
            for s in J: cols[s] = rows[s]
            keep = [s for s in range(n) if s not in J]
            outs = ''.join(rows[s] for s in keep) + ''.join(cols[s] for s in keep)
            R = np.einsum(''.join(rows) + ''.join(cols) + '->' + outs, T)
            tot += mpq(-1, 2) ** k * (sum(x * x for x in np.ravel(R)) if np.ndim(R) else R * R)
    return tot
