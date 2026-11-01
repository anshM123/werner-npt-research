"""Four-replica SDP relaxation of the projector case of the endpoint (d = 3):
   min over rank-RANK projectors P on (C^3)^(x)n of q_n(P) = Tr[(P (x) P) G],  G = (x)_i (F_i - I/2).
Variable Z = P^(x)4 twirled by local U(3)^n;  Z lies in A^(x)n,  A = span{sigma : sigma in S_4} on (C^3)^(x)4 per site.
Z = (+)_patterns Z_lam (x) I, with Z_lam real symmetric on (x)_s S^{lam_s}; we use hat blocks Zhat_lam = (prod m) Z_lam.
Symmetry: Z commutes with the simultaneous replica action of S_4 and is invariant under site permutations, so for each
sorted pattern lam, Zhat_lam = sum_c U_c (W_c (x) I_{d_c}) U_c^T with W_c symmetric (isotypic decomposition of
S_4 x Stab(lam) on (x)_s S^{lam_s}).  Z >= 0  <=>  all W_c >= 0.

Constraints (all valid for twirl(P^(x)4); replicas 0..3):
   moments  Tr[Z tau^glob] = RANK^{cycles(tau)}                                         (always)
   Q<SQ>T<ST> : twirl( (x)_r A_r^{(T if r in ST)} ) >= 0 with A_r = I - P for r in SQ, A_r = P otherwise.
                Expanding I - P gives sum_{T subset SQ} (-1)^{|SQ|-|T|} [Tr_T Z (x) I_T] / RANK^{|T|}, then the
                partial transpose on ST.  ST must be '' (target A), '3' (target B_{3,1}) or '23' (target B_{2,2}).
   aliases:  LE = Q3T (P <= I),  PT1 = QT3,  PT2 = QT23,  PT1LE = Q3T3,  PT2LE = Q3T23,
             QQ = Q23T, QQQ = Q123T, QQQQ = Q0123T
   PROJ  : Tr_3[Z F_23^glob] = Tr_3[Z] / RANK   (operator identity on replicas 0,1,2, from P^2 = P)
Objective Tr[Z (G (x) I (x) I)] = RANK^2 q_n(P).
"""
import sys, time, itertools, argparse, numpy as np, scipy.sparse as sp
from alg4 import *

LET_I = 'abcdefgh'; LET_J = 'ijklmnop'; LET_K = 'ABCDEFGH'; LET_L = 'IJKLMNOP'
ALIAS = {'LE': 'Q3T', 'PT1': 'QT3', 'PT2': 'QT23', 'PT1LE': 'Q3T3', 'PT2LE': 'Q3T23',
         'QQ': 'Q23T', 'QQQ': 'Q123T', 'QQQQ': 'Q0123T'}

def kron_list(mats):
    M = np.array([[1.0]])
    for X in mats: M = np.kron(M, X)
    return M

def orbit_of(lam):
    """distinct site permutations of a sorted pattern: list of (lam', perm) with lam'[q] = lam[perm[q]]"""
    n = len(lam); seen = {}
    for perm in itertools.permutations(range(n)):
        lp = tuple(lam[perm[q]] for q in range(n))
        if lp not in seen: seen[lp] = perm
    return list(seen.items())

def stabilizer(lam):
    n = len(lam)
    return [p for p in itertools.permutations(range(n)) if all(lam[p[q]] == lam[q] for q in range(n))]

def svec_pairs(m):
    return [(i, j) for j in range(m) for i in range(j + 1)]

def svec(V):
    m = V.shape[-1]
    out = [V[..., i, j] * (1.0 if i == j else np.sqrt(2.0)) for (i, j) in svec_pairs(m)]
    return np.stack(out, axis=-1)

def unsvec(v, m):
    V = np.zeros((m, m)); k = 0
    for j in range(m):
        for i in range(j + 1):
            if i == j: V[i, i] = v[k]
            else: V[i, j] = V[j, i] = v[k] / np.sqrt(2.0)
            k += 1
    return V

def sym_adapt(rep_mats, stab, shape, seed=11, tol=1e-8, tries=6):
    """isotypic decomposition of (replica group, given as the full list rep_mats) x (site permutations stab) acting on
    R^shape.  Returns classes [dict(d, m, U)], U: D x (m d) with columns (copy a, vector i), such that every commutant
    element X satisfies U_c^T X U_c = W_c (x) I_d and U_c^T X U_c' = 0 (c != c')."""
    n = len(shape); D = int(np.prod(shape))
    if D == 1:
        return [dict(d=1, m=1, U=np.ones((1, 1)))]
    def reyn(X):
        Y = sum(R @ X @ R.T for R in rep_mats) / len(rep_mats) if rep_mats else X
        if len(stab) > 1:
            T = Y.reshape(shape + shape); acc = np.zeros_like(T)
            for p in stab: acc += T.transpose(list(p) + [n + a for a in p])
            Y = (acc / len(stab)).reshape(D, D)
        return Y
    for attempt in range(tries):
        rng = np.random.default_rng(seed + 7919 * attempt)
        X = reyn(rng.standard_normal((D, D))); X = (X + X.T) / 2
        ev, V = np.linalg.eigh(X)
        sc = np.abs(ev).max()
        spaces = []; st = 0
        for k in range(1, D + 1):
            if k == D or ev[k] - ev[k - 1] > tol * sc:
                spaces.append(V[:, st:k]); st = k
        Y = reyn(rng.standard_normal((D, D))); ny = np.linalg.norm(Y)
        ns = len(spaces); par = list(range(ns))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for i in range(ns):
            for j in range(i + 1, ns):
                if spaces[i].shape[1] != spaces[j].shape[1]: continue
                if np.linalg.norm(spaces[i].T @ Y @ spaces[j]) > 1e-7 * ny:
                    par[find(i)] = find(j)
        groups = {}
        for i in range(ns): groups.setdefault(find(i), []).append(i)
        classes = []; ok = True
        for g in groups.values():
            dd = spaces[g[0]].shape[1]
            B1 = spaces[g[0]]; cols = [B1]
            for i in g[1:]:
                T = spaces[i].T @ Y @ B1
                s2 = np.trace(T.T @ T) / dd
                if s2 < (1e-6 * ny) ** 2 or not np.allclose(T.T @ T, s2 * np.eye(dd), atol=1e-7 * s2):
                    ok = False; break
                cols.append(spaces[i] @ T / np.sqrt(s2))
            if not ok: break
            classes.append(dict(d=dd, m=len(g), U=np.hstack(cols)))
        if not ok: continue
        Uall = np.hstack([c['U'] for c in classes])
        if not np.allclose(Uall.T @ Uall, np.eye(D), atol=1e-9): continue
        Xv = reyn(rng.standard_normal((D, D))); Xv = Xv + Xv.T
        Yv = Uall.T @ Xv @ Uall; rec = np.zeros_like(Yv); off = 0
        for c in classes:
            k = c['m'] * c['d']
            blk = Yv[off:off + k, off:off + k].reshape(c['m'], c['d'], c['m'], c['d'])
            W = np.einsum('aibi->ab', blk) / c['d']
            rec[off:off + k, off:off + k] = np.kron(W, np.eye(c['d'])); off += k
        err = np.abs(Yv - rec).max() / max(1.0, np.abs(Yv).max())
        if err < 1e-9:
            return classes
    raise RuntimeError(f"sym_adapt failed for shape {shape}")

def reduce_blocks(Img, c):
    """Img (K, D, D) commutant elements -> reduced blocks (K, m, m)"""
    X = np.einsum('dp,kde,eq->kpq', c['U'], Img, c['U'], optimize=True)
    X = X.reshape(Img.shape[0], c['m'], c['d'], c['m'], c['d'])
    return np.einsum('kaibi->kab', X) / c['d']

def contract_sites(Zt, Ts, n):
    sub = 'z' + LET_I[:n] + LET_J[:n] + ',' + ','.join(LET_K[s] + LET_L[s] + LET_I[s] + LET_J[s] for s in range(n)) \
          + '->z' + LET_K[:n] + LET_L[:n]
    return np.einsum(sub, Zt, *Ts, optimize='greedy')

def tr_embed(X, reps):
    """per-site operator on 4 replicas: trace the replicas in reps and put identities back on them"""
    reps = sorted(set(reps))
    if not reps: return X
    if len(reps) == 4: return np.trace(X) * np.eye(DIM)
    T = X.reshape((d,) * 8)
    row = [chr(97 + r) for r in range(4)]; col = [chr(101 + r) for r in range(4)]
    colt = list(col)
    for r in reps: colt[r] = row[r]
    keep = [r for r in range(4) if r not in reps]
    out = ''.join(row[r] for r in keep) + ''.join(col[r] for r in keep)
    Y = np.einsum(''.join(row) + ''.join(colt) + '->' + out, T)
    subs = [out] + [row[r] + col[r] for r in reps]
    ops = [Y] + [np.eye(d)] * len(reps)
    return np.einsum(','.join(subs) + '->' + ''.join(row) + ''.join(col), *ops).reshape(DIM, DIM)

def cycle_lengths(s):
    seen = set(); out = []
    for i in range(len(s)):
        if i not in seen:
            L = 0; j = i
            while j not in seen: seen.add(j); j = s[j]; L += 1
            out.append(L)
    return out

def parse_spec(name):
    """Q<SQ>T<ST>[R1|R2]"""
    name = ALIAS.get(name, name)
    R = 0
    if name.endswith('R1'): R = 1; name = name[:-2]
    elif name.endswith('R2'): R = 2; name = name[:-2]
    assert name.startswith('Q') and 'T' in name, name
    sq, st = name[1:].split('T')
    return tuple(sorted(int(c) for c in sq)), tuple(sorted(int(c) for c in st)), R

def slot_perm(X, order):
    """reorder the four replica slots of a per-site operator: new slot k = old slot order[k]"""
    return X.reshape((d,) * 8).transpose(list(order) + [4 + o for o in order]).reshape(DIM, DIM)

def realign1(X):
    """covariance of the pair (2 ket, 3 bra): Y[(a0,a1,a2,b2),(b0,b1,b3,a3)] = X[(a0,a1,a2,a3),(b0,b1,b2,b3)]; barred slot 3"""
    return X.reshape((d,) * 8).transpose(0, 1, 2, 6, 4, 5, 7, 3).reshape(DIM, DIM)

def realign2(X):
    """covariance of the pairs (0 ket, 1 bra), (2 ket, 3 bra), slots ordered (a0, a2, b0, b2 | b1, b3, a1, a3); barred {2,3}"""
    return X.reshape((d,) * 8).transpose(0, 2, 4, 6, 5, 7, 1, 3).reshape(DIM, DIM)

def constraint_terms(name, RANK, ratio=None):
    """returns (barred-slot key of the target algebra, [(coef, key, per-site dense map)])"""
    SQ, ST, R = parse_spec(name)
    if R == 0:
        barred = set(ST)
    elif R == 1:
        assert not set(ST) & {2, 3} and ((2 in SQ) == (3 in SQ)), name
        barred = set(ST) | {3}
    else:
        assert not ST and ((0 in SQ) == (1 in SQ)) and ((2 in SQ) == (3 in SQ)), name
        barred = {2, 3}
    assert 1 <= len(barred) <= 2 or len(barred) == 0, name
    order = [r for r in range(4) if r not in barred] + sorted(barred)
    tgkey = tuple(range(4 - len(barred), 4))
    def post(Y):
        if R == 1: Y = realign1(Y)
        elif R == 2: Y = realign2(Y)
        if order != [0, 1, 2, 3]: Y = slot_perm(Y, order)
        return Y
    terms = []
    for k in range(len(SQ) + 1):
        for T in itertools.combinations(SQ, k):
            key = f"tr{''.join(map(str, T))}_pt{''.join(map(str, ST))}_R{R}_o{''.join(map(str, order))}"
            phi = (lambda X, T=T: post(pt(tr_embed(X, T), list(ST))))
            r = (1.0 / float(RANK)) if ratio is None else ratio
            terms.append(((-1.0) ** (len(SQ) - k) * r ** k, key, phi))
    return tgkey, terms

def tmap_matrix(Tdict, src, tgt):
    os_ = np.cumsum([0] + [f * f for f in src.f]); ot = np.cumsum([0] + [f * f for f in tgt.f])
    M = np.zeros((ot[-1], os_[-1]))
    for (mu, lam), T in Tdict.items():
        g, f = tgt.f[mu], src.f[lam]
        M[ot[mu]:ot[mu] + g * g, os_[lam]:os_[lam] + f * f] = T.reshape(g * g, f * f)
    return M

def conj_matrix(rhos):
    import scipy.linalg
    return scipy.linalg.block_diag(*[np.kron(r, r) for r in rhos])

def detect_group(Ms, src, tgt, tgkey):
    """slot permutations tau (preserving the barred set) such that for some replica permutation sigma every term map
    intertwines conj(sigma) on the source with conj(tau) on the target"""
    barred = set(tgkey)
    CA = {s: conj_matrix(src.coeffs(PERM[s])) for s in S4}
    grp = []
    for tau in S4:
        if set(tau[r] for r in barred) != barred: continue
        Ct = conj_matrix(tgt.coeffs(PERM[tau]))
        for s in S4:
            if all(np.abs(M @ CA[s] - Ct @ M).max() < 1e-9 for M in Ms):
                grp.append(tau); break
    return grp

class Relax:
    def __init__(self, n, RANK, algs, cons, verbose=True, w=0.5, spec=None):
        self.n = n; self.RANK = RANK; self.cons = cons; self.verbose = verbose; self.w = w
        # spectrum of the operator on its support: projector of rank RANK by default
        self.spec = [1.0] * int(round(RANK)) if spec is None else [float(x) for x in spec]
        pk = lambda k: float(sum(x ** k for x in self.spec))
        self.TR = pk(1); self.MU = max(self.spec)
        A, B1, B2, A3 = algs
        self.A, self.B1, self.B2, self.A3 = algs
        t0 = time.time()
        self.rhoA = {tau: A.coeffs(PERM[tau]) for tau in S4}
        # ---------------- variables ----------------
        self.src = []; off = 0; self.wcones = []
        for lam in itertools.combinations_with_replacement(range(A.nb), n):
            shape = tuple(A.f[l] for l in lam); D = int(np.prod(shape))
            rep_mats = [kron_list([self.rhoA[tau][l] for l in lam]) for tau in S4]
            classes = sym_adapt(rep_mats, stabilizer(lam), shape)
            basis = []
            for c in classes:
                Ua = c['U'].reshape(D, c['m'], c['d'])
                self.wcones.append((off + len(basis), c['m']))
                for (i, j) in svec_pairs(c['m']):
                    Bij = Ua[:, i, :] @ Ua[:, j, :].T
                    basis.append(Bij if i == j else (Bij + Bij.T) / np.sqrt(2.0))
            Bt = np.array(basis)
            self.src.append(dict(lam=lam, shape=shape, D=D, classes=classes, orbit=orbit_of(lam), off=off,
                                 K=len(basis), Bt=Bt))
            off += len(basis)
        self.N = off
        if verbose:
            print(f"[n={n} rank={RANK:g}] {len(self.src)} sorted Z-patterns, {self.N} free parameters, "
                  f"{len(self.wcones)} W-blocks (max size {max(m for _, m in self.wcones)}); {time.time() - t0:.1f}s", flush=True)
        # ---------------- linear functionals ----------------
        self.obj = self.functional(PERM[(1, 0, 2, 3)] - w * np.eye(DIM))
        self.mom_rows = []; self.mom_rhs = []
        for tau in [(0, 1, 2, 3), (1, 0, 2, 3), (1, 2, 0, 3), (1, 2, 3, 0), (1, 0, 3, 2)]:
            self.mom_rows.append(self.functional(PERM[tau])); self.mom_rhs.append(float(np.prod([pk(L) for L in cycle_lengths(tau)])))
        # ---------------- constraint blocks ----------------
        self.maps = {}
        self.blocks = []   # list of (name, mu, class, M (svec rows x N))
        algmap = {(): A, (3,): B1, (2, 3): B2}
        for cname in cons:
            if cname == 'PROJ': continue
            if cname == 'CUB': continue
            tgkey, tdefs = constraint_terms(cname, RANK, ratio=self.MU / self.TR)
            tg = algmap[tgkey]
            terms = []
            for coef, key, phi in tdefs:
                if key not in self.maps: self.maps[key] = site_map(phi, A, tg)
                terms.append((coef, key))
            grp = detect_group([tmap_matrix(self.maps[k], A, tg) for _, k in terms], A, tg, tgkey)
            nb0 = len(self.blocks)
            self.add_psd(cname, tg, grp, terms)
            if verbose:
                nb = self.blocks[nb0:]
                print(f"   {cname} ({len(terms)} terms, target {tg.name}, symmetry group {len(grp)}): {len(nb)} reduced blocks, "
                      f"max size {max(b[2]['m'] for b in nb)}, rows {sum(b[3].shape[0] for b in nb)}; {time.time() - t0:.1f}s", flush=True)
        self.eq_rows = np.zeros((0, self.N)); self.eq_rhs = np.zeros(0)
        if 'PROJ' in cons:
            F23 = PERM[(0, 1, 3, 2)]
            self.maps['pa'] = site_map(lambda X: ptrace_last(X @ F23), A, A3)
            self.maps['pb'] = site_map(lambda X: ptrace_last(X), A, A3)
            assert max(self.spec) - min(self.spec) < 1e-12, 'PROJ needs a multiple of a projector'
            rows = self.linear_image(A3, [(1.0, 'pa'), (-self.MU / self.TR, 'pb')])
            u, s, vt = np.linalg.svd(rows, full_matrices=False)
            rk = int(np.sum(s > 1e-10 * s[0]))
            self.eq_rows = vt[:rk]; self.eq_rhs = np.zeros(rk)
            if verbose: print(f"   PROJ: {rows.shape[0]} identity rows, rank {rk}; {time.time() - t0:.1f}s", flush=True)
        if 'CUB' in cons:
            # minimal polynomial x (x - l1)(x - l2) of a rank-2 operator with eigenvalues l1, l2 on its support:
            # Tr_23[Z pi_123] - (e1/TR) Tr_23[Z F_12] + (e2/TR^2) Tr_23[Z] = 0  (operator identity on replicas 0,1)
            ev = sorted(set(round(x, 14) for x in self.spec)); assert len(ev) <= 2
            e1 = sum(ev); e2 = ev[0] * ev[1] if len(ev) == 2 else 0.0
            if len(ev) == 1: e1, e2 = 2 * ev[0], ev[0] ** 2      # (x - l)^2 x, implied by PROJ in that case
            A2 = Alg([perm_op(s, 2) for s in [(0, 1), (1, 0)]], 'A2')
            pi = PERM[(0, 2, 3, 1)]; F12 = PERM[(0, 2, 1, 3)]
            tr23 = lambda X: ptrace_last(ptrace_last(X), nr=3)
            self.maps['cub_a'] = site_map(lambda X: tr23(X @ pi), A, A2)
            self.maps['cub_b'] = site_map(lambda X: tr23(X @ F12), A, A2)
            self.maps['cub_c'] = site_map(lambda X: tr23(X), A, A2)
            rows = self.linear_image(A2, [(1.0, 'cub_a'), (-e1 / self.TR, 'cub_b'), (e2 / self.TR ** 2, 'cub_c')])
            u, s, vt = np.linalg.svd(rows, full_matrices=False)
            rk = int(np.sum(s > 1e-10 * s[0]))
            self.eq_rows = np.vstack([self.eq_rows, vt[:rk]]); self.eq_rhs = np.concatenate([self.eq_rhs, np.zeros(rk)])
            if verbose: print(f"   CUB: {rows.shape[0]} identity rows, rank {rk}; {time.time() - t0:.1f}s", flush=True)
        self.build_time = time.time() - t0

    # ------------------------------------------------------------------
    def functional(self, X):
        """coefficients of w -> Tr[Z (X^(x)n)] for a per-site element X of A"""
        c1 = site_functional(X, self.A)
        out = np.zeros(self.N)
        for s in self.src:
            C = kron_list([c1[l] for l in s['lam']])
            out[s['off']:s['off'] + s['K']] = len(s['orbit']) * np.einsum('kIJ,IJ->k', s['Bt'], C)
        return out

    def images(self, mu, terms, Dmu):
        """for target pattern mu: list over sorted source patterns of Img (K_lam, Dmu, Dmu)"""
        n = self.n; out = []
        for s in self.src:
            K = s['K']; Img = np.zeros((K, Dmu, Dmu)); nz = False
            Zfull = s['Bt'].reshape((K,) + s['shape'] + s['shape'])
            for lp, perm in s['orbit']:
                Zt = None
                for coef, mname in terms:
                    Tm = self.maps[mname]
                    Ts = [Tm[(mu[q], lp[q])] for q in range(n)]
                    if any(not np.any(T) for T in Ts): continue
                    if Zt is None:
                        Zt = Zfull.transpose([0] + [1 + p for p in perm] + [1 + n + p for p in perm])
                    Img += coef * contract_sites(Zt, Ts, n).reshape(K, Dmu, Dmu); nz = True
            out.append(Img if nz else None)
        return out

    def add_psd(self, cname, tg, grp, terms):
        n = self.n
        rho = {tau: tg.coeffs(PERM[tau]) for tau in grp}
        for mu in itertools.combinations_with_replacement(range(tg.nb), n):
            shape = tuple(tg.f[l] for l in mu); Dmu = int(np.prod(shape))
            rep_mats = [kron_list([rho[tau][l] for l in mu]) for tau in grp]
            classes = sym_adapt(rep_mats, stabilizer(mu), shape)
            imgs = self.images(mu, terms, Dmu)
            for c in classes:
                M = np.zeros((c['m'] * (c['m'] + 1) // 2, self.N))
                for s, Img in zip(self.src, imgs):
                    if Img is None: continue
                    V = reduce_blocks(Img, c)
                    V = (V + V.transpose(0, 2, 1)) / 2
                    M[:, s['off']:s['off'] + s['K']] = svec(V).T
                if np.abs(M).max() > 1e-12:
                    self.blocks.append((cname, mu, c, M))

    def linear_image(self, tg, terms):
        rows = []
        for mu in itertools.combinations_with_replacement(range(tg.nb), self.n):
            Dmu = int(np.prod([tg.f[l] for l in mu]))
            imgs = self.images(mu, terms, Dmu)
            R = np.zeros((Dmu * Dmu, self.N))
            for s, Img in zip(self.src, imgs):
                if Img is not None: R[:, s['off']:s['off'] + s['K']] = Img.reshape(s['K'], -1).T
            rows.append(R)
        return np.vstack(rows)

    # ------------------------------------------------------------------
    def solve(self, tol=1e-9, max_iter=500, verbose=False, extra_eq=None):
        import clarabel
        N = self.N
        Arows = []; b = []; cones = []
        E = [np.array(self.mom_rows), self.eq_rows]; Eb = list(self.mom_rhs) + list(self.eq_rhs)
        if extra_eq is not None:
            E.append(extra_eq[0]); Eb += list(extra_eq[1])
        Arows.append(np.vstack(E)); b += Eb
        cones.append(clarabel.ZeroConeT(len(Eb)))
        nn_rows = []
        for (o, m) in self.wcones:
            if m == 1:
                r = np.zeros(N); r[o] = -1.0; nn_rows.append(r)
        for (cname, mu, c, M) in self.blocks:
            if c['m'] == 1: nn_rows.append(-M[0])
        if nn_rows:
            Arows.append(np.array(nn_rows)); b += [0.0] * len(nn_rows); cones.append(clarabel.NonnegativeConeT(len(nn_rows)))
        for (o, m) in self.wcones:
            if m > 1:
                k = m * (m + 1) // 2; Rr = np.zeros((k, N)); Rr[np.arange(k), o + np.arange(k)] = -1.0
                Arows.append(Rr); b += [0.0] * k; cones.append(clarabel.PSDTriangleConeT(m))
        for (cname, mu, c, M) in self.blocks:
            if c['m'] > 1:
                Arows.append(-M); b += [0.0] * M.shape[0]; cones.append(clarabel.PSDTriangleConeT(c['m']))
        Amat = sp.csc_matrix(np.vstack(Arows)); bvec = np.array(b)
        P = sp.csc_matrix((N, N))
        st = clarabel.DefaultSettings(); st.verbose = verbose; st.max_iter = max_iter
        st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
        t0 = time.time()
        solver = clarabel.DefaultSolver(P, self.obj.copy(), Amat, bvec, cones, st)
        sol = solver.solve()
        self.sol = sol; self.Amat = Amat; self.bvec = bvec; self.cones = cones
        self.solve_time = time.time() - t0
        return sol

    # ------------------------------------------------------------------
    def evaluate(self, w):
        """min eigenvalue per constraint family (W blocks = 'Z'), moment and equality residuals, objective"""
        mins = {}
        for (o, m) in self.wcones:
            k = m * (m + 1) // 2
            W = unsvec(w[o:o + k], m); mins['Z'] = min(mins.get('Z', np.inf), np.linalg.eigvalsh(W).min())
        for (cname, mu, c, M) in self.blocks:
            V = unsvec(M @ w, c['m']); mins[cname] = min(mins.get(cname, np.inf), np.linalg.eigvalsh(V).min())
        mom = np.array(self.mom_rows) @ w - np.array(self.mom_rhs)
        eq = self.eq_rows @ w if self.eq_rows.shape[0] else np.zeros(0)
        return mins, np.abs(mom).max(), (np.abs(eq).max() if eq.size else 0.0), float(self.obj @ w)

    def w_from_zfull(self, Zfull):
        """project a full hat-coefficient tensor (23,)*n onto the parameters; returns w and the max residual of
        the reconstruction (zero iff Zfull is invariant and symmetric)"""
        A = self.A; offs = np.cumsum([0] + [f * f for f in A.f])
        w = np.zeros(self.N); res = 0.0
        for s in self.src:
            lam = s['lam']; n = len(lam)
            Zl = Zfull[np.ix_(*[np.arange(offs[l], offs[l + 1]) for l in lam])]
            Zl = Zl.reshape([A.f[l] for l in lam for _ in (0, 1)])            # (i0 j0 i1 j1 ...)
            Zl = Zl.transpose([2 * q for q in range(n)] + [2 * q + 1 for q in range(n)]).reshape(s['D'], s['D'])
            rec = np.zeros_like(Zl); k0 = 0
            for c in s['classes']:
                V = reduce_blocks(Zl[None], c)[0]
                rec += c['U'] @ np.kron(V, np.eye(c['d'])) @ c['U'].T
                vals = svec((V + V.T) / 2)
                w[s['off'] + k0: s['off'] + k0 + len(vals)] = vals; k0 += len(vals)
            res = max(res, np.abs(rec - Zl).max())
        return w, res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int); ap.add_argument('--rank', type=float, default=2.0)
    ap.add_argument('--cons', default='LE,PT1,PT2')
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('--save', default='')
    ap.add_argument('--w', type=float, default=0.5)
    ap.add_argument('--spec', default='')
    a = ap.parse_args()
    algs = build_algebras()
    cons = [c for c in a.cons.split(',') if c]
    spec = [float(x) for x in a.spec.split(',')] if a.spec else None
    rel = Relax(a.n, a.rank, algs, cons, w=a.w, spec=spec)
    sol = rel.solve(tol=a.tol, verbose=a.verbose)
    R = rel.TR
    print(f"n={a.n} spec={rel.spec} w={a.w:g} cons={'+'.join(cons)}: status {sol.status}, iters {sol.iterations}, "
          f"obj {sol.obj_val:.10e} -> q_n relaxation = {sol.obj_val / R ** 2:.10e}  (build {rel.build_time:.1f}s, solve {rel.solve_time:.1f}s)", flush=True)
    w = np.array(sol.x)
    mins, mres, eres, objv = rel.evaluate(w)
    print("   min eig per constraint:", {k: f"{v:.2e}" for k, v in mins.items()}, f" moment res {mres:.1e}, eq res {eres:.1e}")
    if a.save:
        np.savez(a.save, w=w, z=np.array(sol.z), s=np.array(sol.s), obj=sol.obj_val, n=a.n, rank=R, cons=a.cons)
