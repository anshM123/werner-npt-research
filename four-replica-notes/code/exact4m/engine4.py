"""Exact linear data of the two-letter relaxation for n sites via a shared tree of per-site maps.

Per-site maps (exact_core.SiteMap) are scaled to integer matrices T[(mu_s, lam_s)] (g^2 x f^2, common denominator
dT).  For a source group (c, lam) with integer tensor Xt (K, D, D) in pair layout (K, f_0^2, ..., f_{n-1}^2), the
images in ALL unsorted target patterns mu' are produced by a tree over the sites (prefixes shared); the image of the
site orbit of lam in a sorted target mu is the sum over orbit representatives (lam', perm) of the image of lam in
mu' = mu o perm^{-1}, with the site axes permuted back.  All tree arithmetic is in float64 on integers with an
a-priori bound < 2^53 (checked), hence exact.
"""
import itertools, math, numpy as np, flint
from core4 import (site_algebras, FAMS, letter_group, locate, exact_coeffs, SiteMap, fm_den, fm_to_int, fm_kron,
                   perm_matrix_fm, sorted_patterns, orbit_perms, stabilizer, S4, PERM, ptrace_last, inv, realign1, fq)

EXACT_LIMIT = 2 ** 53

# ----------------------------------------------------------------------------- variable space for n sites
class VarSpace:
    def __init__(self, A, n, log=print):
        import time
        t0 = time.time(); self.A = A; self.n = n; nb = len(A.blocks)
        self.f = [b['f'] for b in A.blocks]; self.m = [b['m'] for b in A.blocks]
        self.N = 0; self.groups = {}; self.order = []
        for c in FAMS:
            Gc = letter_group(c)
            for lam in sorted_patterns(nb, n):
                shape = tuple(self.f[l] for l in lam); D = int(np.prod(shape))
                Pi = flint.fmpq_mat(D, D)
                for s in Gc: Pi = Pi + fm_kron([A.blocks[l]['rho'][s] for l in lam])
                Pi = Pi * fq(1, len(Gc))
                R, t = Pi.rref()
                if t == 0: continue
                V = flint.fmpq_mat(D, t, [R[k, j] for j in range(D) for k in range(t)])
                stab = [p for p in stabilizer(lam) if list(p) != list(range(n))]
                VtVi = (V.transpose() * V).inv()
                Ks = []
                for p in stab:
                    Pq = perm_matrix_fm(shape, p)
                    K = VtVi * (V.transpose() * (Pq * V))
                    assert Pq * V == V * K, "site permutation does not preserve range(Pi^T)"
                    Ks.append(K)
                rows = []
                for a in range(t):
                    for b_ in range(a, t):
                        E = flint.fmpq_mat(t, t); E[a, b_] = 1; E[b_, a] = 1
                        Y = E
                        for K in Ks: Y = Y + K * E * K.transpose()
                        rows.append([Y[i, j] for i in range(t) for j in range(t)])
                RR, kY = flint.fmpq_mat(len(rows), t * t, [v for r in rows for v in r]).rref()
                Ys = []; piv = []
                for k in range(kY):
                    Y = flint.fmpq_mat(t, t, [RR[k, e] for e in range(t * t)])
                    e0 = next(e for e in range(t * t) if RR[k, e] != 0)
                    Ys.append(Y); piv.append(divmod(e0, t))
                idx = list(range(self.N, self.N + kY)); self.N += kY
                Gl = fm_kron([A.blocks[l]['G'] for l in lam]); Gil = fm_kron([A.blocks[l]['Gi'] for l in lam])
                Xh = [Gil * (V * Y * V.transpose()) for Y in Ys]
                dX = 1
                for M in Xh:
                    dd = fm_den(M); dX = dX * dd // math.gcd(dX, dd)
                Xt = np.array([fm_to_int(M, dX) for M in Xh], dtype=object)
                self.groups[(c, lam)] = dict(c=c, lam=lam, V=V, Ys=Ys, piv=piv, idx=idx, t=t, D=D, shape=shape,
                                             G=Gl, Gi=Gil, K=kY, orbit=orbit_perms(lam),
                                             mlam=int(np.prod([self.m[l] for l in lam])), dX=dX, Xt=Xt,
                                             Xmax=max(abs(int(v)) for v in Xt.ravel()))
                self.order.append((c, lam))
        self.gid = np.zeros(self.N, dtype=np.int64)          # group index of each parameter
        for gi, key in enumerate(self.order):
            for i in self.groups[key]['idx']: self.gid[i] = gi
        log(f"exact variable space n={n}: {self.N} parameters, {len(self.groups)} groups; {time.time() - t0:.1f}s")

    def pair_float(self, key):
        """float64 pair-layout integer tensor (K, f_0^2, ..., f_{n-1}^2) of a group"""
        g = self.groups[key]; K = g['K']; n = self.n; shape = g['shape']
        Y = g['Xt'].reshape((K,) + shape + shape)
        Y = Y.transpose([0] + [a for q in range(n) for a in (1 + q, 1 + n + q)])
        return np.ascontiguousarray(Y.reshape((K,) + tuple(s * s for s in shape))).astype(np.float64)

    # ---------------------------------------------------------------- exact points
    def x_from_hats(self, hats):
        """hats[c][lam] (sorted lam): exact X_hat as fmpq_mat (D x D); returns (x as list of fmpq, ok)"""
        x = [fq(0)] * self.N; ok = True
        for key in self.order:
            g = self.groups[key]; c, lam = key
            B = g['G'] * hats[c][lam]; V = g['V']
            Vp = (V.transpose() * V).inv() * V.transpose()
            Y = Vp * B * Vp.transpose()
            if not (V * Y * V.transpose() == B): ok = False
            acc = flint.fmpq_mat(g['t'], g['t'])
            for k, (a, b_) in enumerate(g['piv']):
                x[g['idx'][k]] = Y[a, b_]; acc = acc + g['Ys'][k] * Y[a, b_]
            if not (acc == Y): ok = False
        nb = len(self.A.blocks)
        for c in FAMS:
            for lam in sorted_patterns(nb, self.n):
                H = hats[c][lam]
                if (c, lam) not in self.groups and any(H[i, j] != 0 for i in range(H.nrows()) for j in range(H.ncols())): ok = False
        return x, ok


# ----------------------------------------------------------------------------- per-site maps
def site_map_int(A, tg, phi):
    sm = SiteMap(A, tg, phi)
    d = 1
    for T in sm.T.values():
        if T is not None:
            dd = fm_den(T); d = d * dd // math.gcd(d, dd)
    Ti = {}
    for key, T in sm.T.items():
        if T is None: continue
        Ti[key] = np.array(fm_to_int(T, d), dtype=np.int64)
    return Ti, d

def functional_map_int(A, w):
    """per-site 'map' to scalars for the functional X -> Tr[X w]: T[(0, l)] = vec(h_l^T) as (1, f^2)"""
    h = exact_coeffs(A, w)
    d = 1
    for hl in h:
        dd = fm_den(hl); d = d * dd // math.gcd(d, dd)
    Ti = {}
    for l, hl in enumerate(h):
        f = hl.nrows(); H = fm_to_int(hl, d)
        v = np.array([[int(H[j, i]) for i in range(f) for j in range(f)]], dtype=np.int64)   # index i*f + j -> h[j, i]
        if np.any(v): Ti[(0, l)] = v
    return Ti, d

def rowsum_bound(Ti):
    return max(int(np.abs(T).sum(axis=1).max()) for T in Ti.values())

# ----------------------------------------------------------------------------- tree
def tree(Xp, lam, Ti, ntg):
    """Xp float64 (K, F_0, ..., F_{n-1}); returns dict mu' -> float64 (K, G_0, ..., G_{n-1}) (G = g^2)"""
    n = len(lam)
    Tf = {k: v.astype(np.float64) for k, v in Ti.items()}
    level = {(): Xp}
    for s in range(n):
        new = {}
        for pre, Y in level.items():
            for ms in range(ntg):
                T = Tf.get((ms, lam[s]))
                if T is None: continue
                new[pre + (ms,)] = np.tensordot(Y, T, axes=([1], [1]))
        level = new
    return level

def assemble(level, orbit, mu):
    """pair-layout image of the site orbit of lam in the sorted target mu: (K, G(mu_0), ..., G(mu_{n-1})) or None"""
    n = len(mu); acc = None
    for lp, perm in orbit:
        pinv = [0] * n
        for q in range(n): pinv[perm[q]] = q
        mup = tuple(mu[pinv[r]] for r in range(n))
        Z = level.get(mup)
        if Z is None: continue
        Y = Z.transpose([0] + [1 + perm[q] for q in range(n)])
        acc = Y.copy() if acc is None else acc + Y
    return acc

def pair_to_mat(acc, g):
    K = acc.shape[0]; n = len(g)
    Y = acc.reshape((K,) + tuple(x for q in range(n) for x in (g[q], g[q])))
    Y = Y.transpose([0] + [1 + 2 * q for q in range(n)] + [2 + 2 * q for q in range(n)])
    Dm = int(np.prod(g))
    return Y.reshape(K, Dm, Dm)

def check_exact(Y, bound):
    assert bound < EXACT_LIMIT, f"float64 bound {bound:.3e} too large"
    assert np.all(np.rint(Y) == Y), "non-integral float64 result"

# ----------------------------------------------------------------------------- map definitions
def cov_entries():
    """list of (xy, r, s, family c, slot permutation P) for the 12 covariance entries"""
    out = []
    for xy in ['uu', 'uv', 'vv']:
        for r, lr in enumerate('uv'):
            for s_, ls in enumerate('uv'):
                c, s = locate(xy + lr + ls)
                out.append((xy, r, s_, c, s))
    return out

def identity_defs():
    """list of (name, [(coef, family, phi)]) for ORTH and MARG (target A3)"""
    out = []
    for c in FAMS:
        for i, j in itertools.combinations(range(4), 2):
            if c[i] == c[j]: continue
            sw = list(range(4)); sw[i], sw[j] = j, i; Fij = PERM[tuple(sw)]
            order = [r for r in range(4) if r != j] + [j]; mv = PERM[tuple(order.index(r) for r in range(4))]
            out.append((f"ORTH{c}{i}{j}", [(1, c, lambda X, F=Fij, M=mv: ptrace_last(M @ (X @ F) @ M.T))]))
    for xyz in ['uuu', 'uuv', 'uvv', 'vvv']:
        terms = []
        for last, sg in (('u', 1), ('v', -1)):
            c, s = locate(xyz + last); P = PERM[s]
            terms.append((sg, c, lambda X, P=P: ptrace_last(P @ X @ P.T)))
        out.append((f"MARG{xyz}", terms))
    return out
