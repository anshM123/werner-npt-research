"""Two-letter four-replica relaxation for Hermitian rank-2 operators C = cu |u><u| + cv |v><v| (u, v orthonormal) on
(C^3)^(x)n, with the replica-swap identities built into the parametrisation.

Variables: X_c = twirl( (x)_i |c_i><c_i| ) for the canonical letter strings c in {uuuu, uuuv, uuvv, uvvv, vvvv}.
Because the letters are pure, X_c F_ij = X_c whenever c_i = c_j (F_ij = global swap of replicas i, j).  In the Fourier
picture this says that every block X_lam is supported on the subspace of (x)_s S^{lam_s} fixed by the letter stabiliser
G_c (diagonal action).  We parametrise X_lam = V (sum_k U_k (W_k (x) I) U_k^T) V^T with V an orthonormal basis of that
subspace and U_k the isotypic bases of the site stabiliser acting on it; X_c >= 0  <=>  all W_k >= 0.
Constraints (each holds for every orthonormal pair u, v):
   moments  Tr[X_c tau^glob] = 1 if every cycle of tau carries one letter, else 0
   ORTH     Tr_j[X_c F_ij] = 0 when c_i != c_j                         (<u|v> = 0)
   MARG     Tr_3 X_{xyzu} = Tr_3 X_{xyzv}                              (unit traces)
   PT       whole-replica partial transposes, one and two replicas, every letter pattern
   COV      [[R(X_xyuu), R(X_xyuv)], [R(X_xyvu), R(X_xyvv)]] >= 0, R = realignment of the replica pair (2,3)
   LE2      Tr_3 X_{xyzu} (x) I - X_{xyzu} - X_{xyzv} >= 0               (|u><u| + |v><v| <= I)
Objective q_n(C) = cu^2 m_uu + cv^2 m_vv + 2 cu cv m_uv,  m_xy = Tr[X_{xy..} (G (x) I (x) I)].
usage: python herm4.py n cu cv [cons]"""
import sys, time, itertools, numpy as np, scipy.sparse as sp
from alg4 import *
from rep4 import (kron_list, orbit_of, stabilizer, svec_pairs, svec, unsvec, sym_adapt, reduce_blocks, contract_sites,
                  tr_embed, realign1, tmap_matrix, conj_matrix)

FAMS = ['uuuu', 'uuuv', 'uuvv', 'uvvv', 'vvvv']
ALL_CONS = ['ORTH', 'MARG', 'PT', 'COV', 'LE2']

def locate(ell):
    c = ''.join(sorted(ell))
    for s in S4:
        if all(ell[s[i]] == c[i] for i in range(4)): return c, s
    raise ValueError(ell)

def letter_group(c):
    return [s for s in S4 if all(c[s[i]] == c[i] for i in range(4))]

def perm_matrix(shape, p):
    """matrix of the site permutation acting on R^shape (new axis q = old axis p[q])"""
    D = int(np.prod(shape))
    E = np.eye(D).reshape((D,) + tuple(shape)).transpose([0] + [1 + a for a in p]).reshape(D, D)
    return E.T

class Herm4:
    def __init__(self, n, algs, cons, verbose=True, w=0.5, pt_cuts=(1, 2)):
        self.n = n; self.w = w; self.verbose = verbose; self.cons = cons
        A, B1, B2, A3 = algs; self.A, self.B1, self.B2, self.A3 = algs
        self.alg = {(): A, (3,): B1, (2, 3): B2}
        t0 = time.time()
        rhoA = {tau: A.coeffs(PERM[tau]) for tau in S4}
        self.fam = {}; off = 0; self.wcones = []
        for c in FAMS:
            G = letter_group(c); srcs = []; start = off
            for lam in itertools.combinations_with_replacement(range(A.nb), n):
                shape = tuple(A.f[l] for l in lam); D = int(np.prod(shape))
                Pi = sum(kron_list([rhoA[tau][l] for l in lam]) for tau in G) / len(G)
                ev, V = np.linalg.eigh((Pi + Pi.T) / 2); V = V[:, ev > 0.5]; t = V.shape[1]
                if t == 0: continue
                mats = [V.T @ perm_matrix(shape, p) @ V for p in stabilizer(lam)]
                classes = [dict(d=cl['d'], m=cl['m'], U=V @ cl['U']) for cl in sym_adapt(mats, [(0,)], (t,))]
                basis = []
                for cl in classes:
                    Ua = cl['U'].reshape(D, cl['m'], cl['d'])
                    self.wcones.append((off + len(basis), cl['m']))
                    for (i, j) in svec_pairs(cl['m']):
                        Bij = Ua[:, i, :] @ Ua[:, j, :].T
                        basis.append(Bij if i == j else (Bij + Bij.T) / np.sqrt(2.0))
                srcs.append(dict(lam=lam, shape=shape, D=D, classes=classes, orbit=orbit_of(lam), off=off,
                                 K=len(basis), Bt=np.array(basis)))
                off += len(basis)
            self.fam[c] = dict(G=G, srcs=srcs, start=start, stop=off)
        self.N = off
        if verbose:
            counts = {c: self.fam[c]['stop'] - self.fam[c]['start'] for c in FAMS}
            print(f"[herm4 n={n}] {self.N} parameters {counts}, {len(self.wcones)} W-blocks "
                  f"(max {max(m for _, m in self.wcones)}); {time.time() - t0:.1f}s", flush=True)
        self.maps = {}
        # ---------------- functionals ----------------
        G01 = PERM[(1, 0, 2, 3)] - w * np.eye(DIM); G12 = PERM[(0, 2, 1, 3)] - w * np.eye(DIM)
        self.m_uu = self.functional('uuuu', G01); self.m_vv = self.functional('vvvv', G01); self.m_uv = self.functional('uuvv', G12)
        rows = []; rhs = []
        for c in FAMS:
            for tau in S4:
                val = 1.0; seen = set()
                for i in range(4):
                    if i in seen: continue
                    j = i; letters = set()
                    while j not in seen: seen.add(j); letters.add(c[j]); j = tau[j]
                    if len(letters) > 1: val = 0.0
                rows.append(self.functional(c, PERM[tau])); rhs.append(val)
        R = np.array(rows); b = np.array(rhs)
        U_, S_, Vt_ = np.linalg.svd(R, full_matrices=False); rk = int(np.sum(S_ > 1e-10 * S_[0]))
        self.mom_rows = np.diag(S_[:rk]) @ Vt_[:rk]; self.mom_rhs = U_[:, :rk].T @ b
        self.mom_full = (R, b)
        assert np.abs(R @ np.linalg.lstsq(R, b, rcond=None)[0] - b).max() < 1e-9
        # ---------------- identities (rows accumulated into an orthonormal basis) ----------------
        self.Q = np.zeros((0, self.N)); self.id_log = []
        if 'ORTH' in cons:
            for c in FAMS:
                for i, j in itertools.combinations(range(4), 2):
                    if c[i] == c[j]: continue
                    sw = list(range(4)); sw[i], sw[j] = j, i; Fij = PERM[tuple(sw)]
                    order = [r for r in range(4) if r != j] + [j]
                    mv = PERM[tuple(order.index(r) for r in range(4))]
                    mk = self.getmap(f'orth_{c}_{i}{j}', lambda X, F=Fij, M=mv: ptrace_last(M @ (X @ F) @ M.T), A3)
                    self.add_identity(f'ORTH{c}{i}{j}', A3, [(1.0, c, mk)])
        if 'MARG' in cons:
            for xyz in ['uuu', 'uuv', 'uvv', 'vvv']:
                t = []
                for last, sg in (('u', 1.0), ('v', -1.0)):
                    c, s = locate(xyz + last); P = PERM[s]
                    t.append((sg, c, self.getmap(f'marg_{xyz}{last}', lambda X, P=P: ptrace_last(P @ X @ P.T), A3)))
                self.add_identity(f'MARG{xyz}', A3, t)
        if verbose:
            print(f"   identities: rank {self.Q.shape[0]} ({'; '.join(self.id_log)}); {time.time() - t0:.1f}s", flush=True)
        # ---------------- PSD constraints ----------------
        self.blocks = []
        if 'PT' in cons:
            for c in FAMS:
                seen = set()
                for k in pt_cuts:
                    for S in itertools.combinations(range(4), k):
                        key = (k, tuple(sorted(c[i] for i in S)))
                        if key in seen: continue
                        seen.add(key)
                        order = [r for r in range(4) if r not in S] + list(S)
                        mv = PERM[tuple(order.index(r) for r in range(4))]
                        cut = [3] if k == 1 else [2, 3]
                        tgkey = (3,) if k == 1 else (2, 3)
                        mk = self.getmap(f'pt_{c}_{S}', lambda X, M=mv, cut=cut: pt(M @ X @ M.T, cut), self.alg[tgkey])
                        self.add_block(f'PT{k}{c}{S}', self.alg[tgkey], tgkey, [[[(1.0, c, mk)]]])
                if verbose: print(f"      PT {c}: {len(self.blocks)} blocks so far; {time.time() - t0:.1f}s", flush=True)
        if 'COV' in cons:
            for xy in ['uu', 'uv', 'vv']:
                ent = [[None, None], [None, None]]
                for r, lr in enumerate('uv'):
                    for s_, ls in enumerate('uv'):
                        c, s = locate(xy + lr + ls); P = PERM[s]
                        mk = self.getmap(f'cov_{xy}{lr}{ls}', lambda X, P=P: realign1(P @ X @ P.T), B1)
                        ent[r][s_] = [(1.0, c, mk)]
                self.add_block(f'COV{xy}', B1, (3,), ent)
            if verbose: print(f"      COV done: {len(self.blocks)} blocks; {time.time() - t0:.1f}s", flush=True)
        if 'LE2' in cons:
            for xyz in ['uuu', 'uuv', 'uvv', 'vvv']:
                terms = []
                c, s = locate(xyz + 'u'); P = PERM[s]
                terms.append((1.0, c, self.getmap(f'le_tr_{xyz}', lambda X, P=P: tr_embed(P @ X @ P.T, [3]), A)))
                for last in 'uv':
                    c2, s2 = locate(xyz + last); P2 = PERM[s2]
                    terms.append((-1.0, c2, self.getmap(f'le_x_{xyz}{last}', lambda X, P=P2: P @ X @ P.T, A)))
                self.add_block(f'LE2{xyz}', A, (), [[terms]])
        self.nrows = sum(b[2] * (b[2] + 1) // 2 for b in self.blocks)
        self.nnz = sum(sum(Mc.size for Mc in b[3].values()) for b in self.blocks)
        if verbose:
            from collections import Counter
            cnt = Counter(b[0][:3] for b in self.blocks)
            print(f"   PSD blocks: {len(self.blocks)} ({dict(cnt)}), max size {max(b[2] for b in self.blocks)}, "
                  f"{self.nrows} svec rows, {self.nnz / 1e6:.1f}M stored entries; {time.time() - t0:.1f}s", flush=True)
        self.build_time = time.time() - t0

    # ------------------------------------------------------------------
    def getmap(self, key, phi, tg):
        if key not in self.maps: self.maps[key] = (site_map(phi, self.A, tg), tg)
        return key

    def functional(self, c, X):
        c1 = site_functional(X, self.A); out = np.zeros(self.N)
        for s in self.fam[c]['srcs']:
            C = kron_list([c1[l] for l in s['lam']])
            out[s['off']:s['off'] + s['K']] = len(s['orbit']) * np.einsum('kIJ,IJ->k', s['Bt'], C)
        return out

    def images(self, mu, terms, Dmu):
        n = self.n; out = {}
        for coef, c, mk in terms:
            Tm = self.maps[mk][0]
            for si, s in enumerate(self.fam[c]['srcs']):
                K = s['K']; Zfull = s['Bt'].reshape((K,) + s['shape'] + s['shape'])
                for lp, perm in s['orbit']:
                    Ts = [Tm[(mu[q], lp[q])] for q in range(n)]
                    if any(not np.any(T) for T in Ts): continue
                    Zt = Zfull.transpose([0] + [1 + p for p in perm] + [1 + n + p for p in perm])
                    img = coef * contract_sites(Zt, Ts, n).reshape(K, Dmu, Dmu)
                    key = (c, si)
                    out[key] = out[key] + img if key in out else img
        return out

    def add_identity(self, name, tg, terms):
        """rows of a linear identity L(X) = 0 (all entries of every target pattern), merged into the orthonormal row
        basis self.Q; reports the singular-value gap used to separate new directions from round-off"""
        kept = 0; gap_lo = np.inf; gap_hi = 0.0
        for mu in itertools.combinations_with_replacement(range(tg.nb), self.n):
            Dmu = int(np.prod([tg.f[l] for l in mu]))
            imgs = self.images(mu, terms, Dmu)
            if not imgs: continue
            R = np.zeros((Dmu * Dmu, self.N))
            for (c, si), Img in imgs.items():
                s = self.fam[c]['srcs'][si]; R[:, s['off']:s['off'] + s['K']] += Img.reshape(s['K'], -1).T
            scale = np.linalg.norm(R, axis=1).max()
            if scale < 1e-12: continue
            R = R / scale
            if self.Q.shape[0]: R = R - (R @ self.Q.T) @ self.Q
            u, sv, vt = np.linalg.svd(R, full_matrices=False)
            good = sv > 1e-8
            if np.any(good): gap_lo = min(gap_lo, sv[good].min())
            if np.any(~good): gap_hi = max(gap_hi, sv[~good].max())
            if np.any(good):
                newQ = vt[good]
                # re-orthonormalise against Q once more for stability
                if self.Q.shape[0]: newQ = newQ - (newQ @ self.Q.T) @ self.Q
                qq, rr = np.linalg.qr(newQ.T); self.Q = np.vstack([self.Q, qq.T]); kept += int(np.sum(good))
        self.id_log.append(f"{name}: +{kept} (kept sv >= {gap_lo:.1e}, dropped <= {gap_hi:.1e})")

    def detect(self, tg, tgkey, ent):
        barred = set(tgkey); grp = []
        CA = {s: conj_matrix(self.A.coeffs(PERM[s])) for s in S4}
        for tau in S4:
            if set(tau[r] for r in barred) != barred: continue
            Ct = conj_matrix(tg.coeffs(PERM[tau])); ok = True
            for row in ent:
                for terms in row:
                    for coef, c, mk in terms:
                        M = tmap_matrix(self.maps[mk][0], self.A, tg)
                        if not any(np.abs(M @ CA[s] - Ct @ M).max() < 1e-9 for s in self.fam[c]['G']): ok = False
            if ok: grp.append(tau)
        return grp

    def add_block(self, name, tg, tgkey, ent):
        n = self.n; R = len(ent)
        grp = self.detect(tg, tgkey, ent)
        rho = {tau: tg.coeffs(PERM[tau]) for tau in grp}
        for mu in itertools.combinations_with_replacement(range(tg.nb), n):
            shape = tuple(tg.f[l] for l in mu); Dmu = int(np.prod(shape))
            imgs = [[self.images(mu, ent[r][s], Dmu) for s in range(R)] for r in range(R)]
            if not any(imgs[r][s] for r in range(R) for s in range(R)): continue
            classes = sym_adapt([kron_list([rho[tau][l] for l in mu]) for tau in grp], stabilizer(mu), shape)
            for cl in classes:
                m = cl['m']; full = {}
                for r in range(R):
                    for s in range(R):
                        for key, Img in imgs[r][s].items():
                            V = reduce_blocks(Img, cl)
                            if key not in full: full[key] = np.zeros((Img.shape[0], R * m, R * m))
                            full[key][:, r * m:(r + 1) * m, s * m:(s + 1) * m] += V
                Md = {}
                for (c, si), F in full.items():
                    F = (F + F.transpose(0, 2, 1)) / 2
                    sdat = self.fam[c]['srcs'][si]
                    if c not in Md: Md[c] = np.zeros(((R * m) * (R * m + 1) // 2, self.fam[c]['stop'] - self.fam[c]['start']))
                    o = sdat['off'] - self.fam[c]['start']
                    Md[c][:, o:o + sdat['K']] += svec(F).T
                Md = {c: Mc for c, Mc in Md.items() if np.abs(Mc).max() > 1e-12}
                if Md: self.blocks.append((name, mu, R * m, Md))

    def block_apply(self, Md, x):
        return sum(Mc @ x[self.fam[c]['start']:self.fam[c]['stop']] for c, Mc in Md.items())

    # ------------------------------------------------------------------
    def assemble(self):
        """sparse constraint matrix and cones in Clarabel form (A x + s = b, s in K); objective set in solve()"""
        import clarabel
        N = self.N
        rows = []; cols = []; vals = []; b = []; cones = []; r0 = 0
        def put_dense(M, colstart):
            nonlocal r0
            ii, jj = np.nonzero(M)
            rows.append(ii + r0); cols.append(jj + colstart); vals.append(M[ii, jj]); r0 += M.shape[0]
        E = np.vstack([self.mom_rows, self.Q]); put_dense(E, 0); b += list(self.mom_rhs) + [0.0] * self.Q.shape[0]
        cones.append(clarabel.ZeroConeT(E.shape[0]))
        nn = 0
        for (o, m) in self.wcones:
            if m == 1:
                rows.append(np.array([r0])); cols.append(np.array([o])); vals.append(np.array([-1.0])); r0 += 1; nn += 1
        for (name, mu_, m, Md) in self.blocks:
            if m == 1:
                for c, Mc in Md.items():
                    ii, jj = np.nonzero(Mc); rows.append(ii + r0); cols.append(jj + self.fam[c]['start']); vals.append(-Mc[ii, jj])
                r0 += 1; nn += 1
        b += [0.0] * nn
        if nn: cones.append(clarabel.NonnegativeConeT(nn))
        for (o, m) in self.wcones:
            if m > 1:
                k = m * (m + 1) // 2
                rows.append(np.arange(k) + r0); cols.append(o + np.arange(k)); vals.append(-np.ones(k)); r0 += k
                b += [0.0] * k; cones.append(clarabel.PSDTriangleConeT(m))
        for (name, mu_, m, Md) in self.blocks:
            if m > 1:
                k = m * (m + 1) // 2
                for c, Mc in Md.items():
                    ii, jj = np.nonzero(Mc); rows.append(ii + r0); cols.append(jj + self.fam[c]['start']); vals.append(-Mc[ii, jj])
                r0 += k; b += [0.0] * k; cones.append(clarabel.PSDTriangleConeT(m))
        A = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(r0, N))
        return A, np.array(b), cones

    def solve(self, cu, cv, tol=1e-9, verbose=False, max_iter=300):
        import clarabel
        if not hasattr(self, 'Aasm'): self.Aasm = self.assemble()
        A, bv, cones = self.Aasm
        N = self.N; c = cu ** 2 * self.m_uu + cv ** 2 * self.m_vv + 2 * cu * cv * self.m_uv
        st = clarabel.DefaultSettings(); st.verbose = verbose; st.max_iter = max_iter
        st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
        t0 = time.time()
        sol = clarabel.DefaultSolver(sp.csc_matrix((N, N)), c, A, bv, cones, st).solve()
        self.sol = sol; self.solve_time = time.time() - t0; self.c = c
        return sol

    def evaluate(self, x, cu, cv):
        mins = {}
        for (o, m) in self.wcones:
            k = m * (m + 1) // 2
            mins['X'] = min(mins.get('X', np.inf), np.linalg.eigvalsh(unsvec(x[o:o + k], m)).min())
        for (name, mu_, m, Md) in self.blocks:
            key = name[:3]
            mins[key] = min(mins.get(key, np.inf), np.linalg.eigvalsh(unsvec(self.block_apply(Md, x), m)).min())
        R, b = self.mom_full
        mom = np.abs(R @ x - b).max()
        eqr = np.abs(self.Q @ x).max() if self.Q.shape[0] else 0.0
        obj = (cu ** 2 * self.m_uu + cv ** 2 * self.m_vv + 2 * cu * cv * self.m_uv) @ x
        return mins, mom, eqr, obj

    def x_from_zfulls(self, zf):
        A = self.A; offs = np.cumsum([0] + [f * f for f in A.f]); x = np.zeros(self.N); res = 0.0
        for c in FAMS:
            done = set()
            for s in self.fam[c]['srcs']:
                lam = s['lam']; n = len(lam); done.add(lam)
                Zl = zf[c][np.ix_(*[np.arange(offs[l], offs[l + 1]) for l in lam])]
                Zl = Zl.reshape([A.f[l] for l in lam for _ in (0, 1)])
                Zl = Zl.transpose([2 * q for q in range(n)] + [2 * q + 1 for q in range(n)]).reshape(s['D'], s['D'])
                rec = np.zeros_like(Zl); k0 = 0
                for cl in s['classes']:
                    V = reduce_blocks(Zl[None], cl)[0]
                    rec += cl['U'] @ np.kron(V, np.eye(cl['d'])) @ cl['U'].T
                    vals = svec((V + V.T) / 2); x[s['off'] + k0: s['off'] + k0 + len(vals)] = vals; k0 += len(vals)
                res = max(res, np.abs(rec - Zl).max())
            # patterns without variables must carry no weight at a true point
            for lam in itertools.combinations_with_replacement(range(A.nb), self.n):
                if lam in done: continue
                Zl = zf[c][np.ix_(*[np.arange(offs[l], offs[l + 1]) for l in lam])]
                res = max(res, np.abs(Zl).max())
        return x, res


def moments_multi(Us, n):
    """m[t] = Tr[(rho_0 (x) rho_1 (x) rho_2 (x) rho_3) t], rho_r = U_r U_r^dagger (site-major convention)"""
    out = np.zeros((24,) * n, complex)
    PT_flat = np.array([PERM[s].T.ravel() for s in S4])
    letters = [chr(c) for c in list(range(65, 91)) + list(range(97, 123))]
    a = letters[:4]; b = letters[4:8]; R = letters[8:12]
    X = [[letters[12 + r * n + s] for s in range(n)] for r in range(4)]
    Ut = [U.reshape((3,) * n + (U.shape[1],)) for U in Us]
    for head in itertools.product(range(24), repeat=n - 1):
        sig = [S4[kk] for kk in head]
        subs = []; ops = []
        for r in range(4):
            subs.append(''.join(X[r][s] for s in range(n - 1)) + a[r] + R[r]); ops.append(Ut[r])
            subs.append(''.join(X[inv(sig[s])[r]][s] for s in range(n - 1)) + b[r] + R[r]); ops.append(Ut[r].conj())
        M = np.einsum(','.join(subs) + '->' + ''.join(a) + ''.join(b), *ops, optimize='greedy').reshape(81, 81)
        out[head] = PT_flat @ M.ravel()
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]); cu = float(sys.argv[2]); cv = float(sys.argv[3])
    cons = sys.argv[4].split(',') if len(sys.argv) > 4 else ALL_CONS
    algs = build_algebras()
    rel = Herm4(n, algs, cons)
    sol = rel.solve(cu, cv)
    x = np.array(sol.x)
    mins, mom, eqr, obj = rel.evaluate(x, cu, cv)
    print(f"herm4 n={n} C = {cu}|u><u| + {cv}|v><v|, cons {'+'.join(cons)}: status {sol.status}, iters {sol.iterations}, "
          f"q relaxation = {sol.obj_val:.10e} (build {rel.build_time:.1f}s, solve {rel.solve_time:.1f}s)")
    print("   min eig:", {k: f"{v:.1e}" for k, v in mins.items()}, f"moment res {mom:.1e}, identity res {eqr:.1e}")
