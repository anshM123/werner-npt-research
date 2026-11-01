"""Two-letter four-replica relaxation for Hermitian (or PSD) rank-2 operators
   C = mu |u><u| + s nu |v><v|,  <u|v> = 0, |u| = |v| = 1, mu, nu > 0, s = -1 (indefinite) or +1 (PSD), on (C^3)^(x)n.
Variables: X_l = twirl( (x)_i |l_i><l_i| ) for letter strings l in {u,v}^4, stored for the canonical strings
uuuu, uuuv, uuvv, uvvv, vvvv (others are replica conjugates).  Each X_c is invariant under the letter stabiliser G_c
and the site permutations, real symmetric; X_c >= 0 is built into the parametrisation.
Valid constraints (all hold for every orthonormal pair u, v):
   moments   Tr[X_l tau^glob] = 1 if every cycle of tau carries one letter, else 0
   SWAP      X_l F_ij = X_l when l_i = l_j            (pure letters)
   ORTH      Tr_j[X_l F_ij] = 0 when l_i != l_j      (<u|v> = 0)
   MARG      Tr_3 X_{xyzu} = Tr_3 X_{xyzv}            (unit traces)
   PT        whole-replica partial transposes X_l^{T_S} >= 0 (one and two replicas, every letter pattern)
   COV       [[R(X_xyuu), R(X_xyuv)], [R(X_xyvu), R(X_xyvv)]] >= 0, R = covariance (realignment) of replicas 2,3:
             the Gram matrix of (x y vec|u><u|, x y vec|v><v|)
   LE2       Tr_3 X_{xyzu} (x) I - X_{xyzu} - X_{xyzv} >= 0   (|u><u| + |v><v| <= I)
Objective q_n(C) = mu^2 m_uu + nu^2 m_vv + 2 s mu nu m_uv, m_xy = Tr[X_{xy..} (G (x) I (x) I)].
usage: python herm2.py n mu nu s [cons]"""
import sys, time, itertools, numpy as np, scipy.sparse as sp
from alg4 import *
from rep4 import (kron_list, orbit_of, stabilizer, svec_pairs, svec, unsvec, sym_adapt, reduce_blocks, contract_sites,
                  tr_embed, realign1, tmap_matrix, conj_matrix)

FAMS = ['uuuu', 'uuuv', 'uuvv', 'uvvv', 'vvvv']

def locate(ell):
    c = ''.join(sorted(ell))
    for s in S4:
        if all(ell[s[i]] == c[i] for i in range(4)): return c, s
    raise ValueError(ell)

def letter_group(c):
    return [s for s in S4 if all(c[s[i]] == c[i] for i in range(4))]

class Herm2:
    def __init__(self, n, algs, cons, verbose=True, w=0.5):
        self.n = n; self.w = w; self.verbose = verbose
        A, B1, B2, A3 = algs; self.A, self.B1, self.B2, self.A3 = algs
        self.alg = {(): A, (3,): B1, (2, 3): B2, 'A3': A3}
        t0 = time.time()
        rhoA = {tau: A.coeffs(PERM[tau]) for tau in S4}
        self.fam = {}; off = 0; self.wcones = []
        for c in FAMS:
            G = letter_group(c); srcs = []
            for lam in itertools.combinations_with_replacement(range(A.nb), n):
                shape = tuple(A.f[l] for l in lam); D = int(np.prod(shape))
                classes = sym_adapt([kron_list([rhoA[tau][l] for l in lam]) for tau in G], stabilizer(lam), shape)
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
            self.fam[c] = dict(G=G, srcs=srcs)
        self.N = off
        if verbose:
            counts = {c: sum(s['K'] for s in self.fam[c]['srcs']) for c in FAMS}
            print(f"[herm2 n={n}] {self.N} parameters {counts}, {len(self.wcones)} W-blocks "
                  f"(max {max(m for _, m in self.wcones)}); {time.time() - t0:.1f}s", flush=True)
        self.maps = {}
        # ---------------- functionals ----------------
        G01 = PERM[(1, 0, 2, 3)] - w * np.eye(DIM); G12 = PERM[(0, 2, 1, 3)] - w * np.eye(DIM)
        self.m_uu = self.functional('uuuu', G01); self.m_vv = self.functional('vvvv', G01); self.m_uv = self.functional('uuvv', G12)
        rows = []; rhs = []
        for c in FAMS:
            for tau in S4:
                val = 1.0
                seen = set()
                for i in range(4):
                    if i in seen: continue
                    j = i; letters = set()
                    while j not in seen: seen.add(j); letters.add(c[j]); j = tau[j]
                    if len(letters) > 1: val = 0.0
                rows.append(self.functional(c, PERM[tau])); rhs.append(val)
        R = np.array(rows); b = np.array(rhs)
        U_, S_, Vt_ = np.linalg.svd(R, full_matrices=False); rk = int(np.sum(S_ > 1e-10 * S_[0]))
        self.mom_rows = np.diag(S_[:rk]) @ Vt_[:rk]; self.mom_rhs = U_[:, :rk].T @ b
        assert np.abs(R @ np.linalg.lstsq(R, b, rcond=None)[0] - b).max() < 1e-9
        # ---------------- equalities ----------------
        eq = []
        if 'SWAP' in cons or 'ORTH' in cons:
            for c in FAMS:
                for i, j in itertools.combinations(range(4), 2):
                    sw = list(range(4)); sw[i], sw[j] = j, i; Fij = PERM[tuple(sw)]
                    if c[i] == c[j] and 'SWAP' in cons:
                        terms = [(1.0, c, self.getmap(f'right_{i}{j}', lambda X, F=Fij: X @ F, A)),
                                 (-1.0, c, self.getmap('id', lambda X: X, A))]
                        eq.append(self.linear_image(A, terms))
                    if c[i] != c[j] and 'ORTH' in cons:
                        # move replica j to slot 3, then trace it after right multiplication by F_ij
                        order = [r for r in range(4) if r != j] + [j]
                        mv = PERM[tuple(order.index(r) for r in range(4))]
                        terms = [(1.0, c, self.getmap(f'orth_{i}{j}', lambda X, F=Fij, M=mv: ptrace_last(M @ (X @ F) @ M.T), A3))]
                        eq.append(self.linear_image(A3, terms))
        if 'MARG' in cons:
            for xyz in itertools.product('uv', repeat=3):
                xyz = ''.join(xyz)
                if xyz != ''.join(sorted(xyz)): continue
                t = []
                for last, sg in (('u', 1.0), ('v', -1.0)):
                    c, s = locate(xyz + last); P = PERM[s]
                    t.append((sg, c, self.getmap(f'marg_{xyz}{last}', lambda X, P=P: ptrace_last(P @ X @ P.T), A3)))
                eq.append(self.linear_image(A3, t))
        if eq:
            E = np.vstack(eq); U_, S_, Vt_ = np.linalg.svd(E, full_matrices=False)
            rk = int(np.sum(S_ > 1e-10 * S_[0])); self.eq_rows = Vt_[:rk]
            if verbose: print(f"   equalities: {E.shape[0]} rows, rank {rk}; {time.time() - t0:.1f}s", flush=True)
        else:
            self.eq_rows = np.zeros((0, self.N))
        # ---------------- PSD constraints ----------------
        self.blocks = []
        if 'PT' in cons:
            for c in FAMS:
                G = self.fam[c]['G']; seen = set()
                for k in (1, 2):
                    for S in itertools.combinations(range(4), k):
                        key = (k, tuple(sorted(c[i] for i in S)))
                        if key in seen: continue
                        seen.add(key)
                        # conjugate so that the cut replicas sit in the last slots, then partial transpose there
                        order = [r for r in range(4) if r not in S] + list(S)
                        mv = PERM[tuple(order.index(r) for r in range(4))]
                        cut = [3] if k == 1 else [2, 3]
                        tgkey = (3,) if k == 1 else (2, 3)
                        mk = self.getmap(f'pt_{c}_{S}', lambda X, M=mv, cut=cut: pt(M @ X @ M.T, cut), self.alg[tgkey])
                        self.add_block(f'PT{c}{S}', self.alg[tgkey], tgkey, [[[(1.0, c, mk)]]])
        if 'COV' in cons:
            for xy in ['uu', 'uv', 'vv']:
                ent = [[None, None], [None, None]]
                for r, lr in enumerate('uv'):
                    for s_, ls in enumerate('uv'):
                        c, s = locate(xy + lr + ls); P = PERM[s]
                        mk = self.getmap(f'cov_{xy}{lr}{ls}', lambda X, P=P: realign1(P @ X @ P.T), B1)
                        ent[r][s_] = [(1.0, c, mk)]
                self.add_block(f'COV{xy}', B1, (3,), ent)
        if 'LE2' in cons:
            for xyz in ['uuu', 'uuv', 'uvv', 'vvv']:
                terms = []
                c, s = locate(xyz + 'u'); P = PERM[s]
                terms.append((1.0, c, self.getmap(f'le_tr_{xyz}', lambda X, P=P: tr_embed(P @ X @ P.T, [3]), A)))
                for last in 'uv':
                    c2, s2 = locate(xyz + last); P2 = PERM[s2]
                    terms.append((-1.0, c2, self.getmap(f'le_x_{xyz}{last}', lambda X, P=P2: P @ X @ P.T, A)))
                self.add_block(f'LE2{xyz}', A, (), [[terms]])
        if verbose:
            from collections import Counter
            cnt = Counter(b[0][:3] for b in self.blocks)
            print(f"   PSD blocks: {len(self.blocks)} ({dict(cnt)}), max size {max(b[2] for b in self.blocks)}; {time.time() - t0:.1f}s", flush=True)
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
        """dict (family, source index) -> Img (K, Dmu, Dmu) for the sum of terms"""
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

    def linear_image(self, tg, terms):
        rows = []
        for mu in itertools.combinations_with_replacement(range(tg.nb), self.n):
            Dmu = int(np.prod([tg.f[l] for l in mu]))
            imgs = self.images(mu, terms, Dmu)
            R = np.zeros((Dmu * Dmu, self.N))
            for (c, si), Img in imgs.items():
                s = self.fam[c]['srcs'][si]; R[:, s['off']:s['off'] + s['K']] += Img.reshape(s['K'], -1).T
            rows.append(R)
        return np.vstack(rows)

    def detect(self, tg, tgkey, ent):
        """replica symmetries (preserving the barred slots) of every entry of a (block) constraint"""
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
            classes = sym_adapt([kron_list([rho[tau][l] for l in mu]) for tau in grp], stabilizer(mu), shape)
            imgs = [[self.images(mu, ent[r][s], Dmu) for s in range(R)] for r in range(R)]
            for cl in classes:
                m = cl['m']; M = np.zeros(((R * m) * (R * m + 1) // 2, self.N))
                # assemble the reduced R m x R m block column by column over parameters
                full = {}
                for r in range(R):
                    for s in range(R):
                        for key, Img in imgs[r][s].items():
                            V = reduce_blocks(Img, cl)
                            if key not in full: full[key] = np.zeros((Img.shape[0], R * m, R * m))
                            full[key][:, r * m:(r + 1) * m, s * m:(s + 1) * m] += V
                for (c, si), F in full.items():
                    F = (F + F.transpose(0, 2, 1)) / 2
                    sdat = self.fam[c]['srcs'][si]
                    M[:, sdat['off']:sdat['off'] + sdat['K']] += svec(F).T
                if np.abs(M).max() > 1e-12: self.blocks.append((name, mu, R * m, M))

    # ------------------------------------------------------------------
    def solve(self, mu, nu, s, tol=1e-9, verbose=False):
        import clarabel
        N = self.N; c = mu ** 2 * self.m_uu + nu ** 2 * self.m_vv + 2 * s * mu * nu * self.m_uv
        Arows = [np.vstack([self.mom_rows, self.eq_rows])]; b = list(self.mom_rhs) + [0.0] * self.eq_rows.shape[0]
        cones = [clarabel.ZeroConeT(len(b))]
        nn = []
        for (o, m) in self.wcones:
            if m == 1: r = np.zeros(N); r[o] = -1.0; nn.append(r)
        for (name, mu_, m, M) in self.blocks:
            if m == 1: nn.append(-M[0])
        if nn: Arows.append(np.array(nn)); b += [0.0] * len(nn); cones.append(clarabel.NonnegativeConeT(len(nn)))
        for (o, m) in self.wcones:
            if m > 1:
                k = m * (m + 1) // 2; Rr = np.zeros((k, N)); Rr[np.arange(k), o + np.arange(k)] = -1.0
                Arows.append(Rr); b += [0.0] * k; cones.append(clarabel.PSDTriangleConeT(m))
        for (name, mu_, m, M) in self.blocks:
            if m > 1: Arows.append(-M); b += [0.0] * M.shape[0]; cones.append(clarabel.PSDTriangleConeT(m))
        A = sp.csc_matrix(np.vstack(Arows)); bv = np.array(b)
        st = clarabel.DefaultSettings(); st.verbose = verbose; st.max_iter = 500
        st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
        t0 = time.time()
        sol = clarabel.DefaultSolver(sp.csc_matrix((N, N)), c, A, bv, cones, st).solve()
        self.sol = sol; self.solve_time = time.time() - t0; self.c = c
        return sol

    def evaluate(self, x, mu, nu, s):
        mins = {}
        for (o, m) in self.wcones:
            k = m * (m + 1) // 2
            mins['X'] = min(mins.get('X', np.inf), np.linalg.eigvalsh(unsvec(x[o:o + k], m)).min())
        for (name, mu_, m, M) in self.blocks:
            key = name[:3]
            mins[key] = min(mins.get(key, np.inf), np.linalg.eigvalsh(unsvec(M @ x, m)).min())
        mom = np.abs(self.mom_rows @ x - self.mom_rhs).max()
        eqr = np.abs(self.eq_rows @ x).max() if self.eq_rows.shape[0] else 0.0
        obj = (mu ** 2 * self.m_uu + nu ** 2 * self.m_vv + 2 * s * mu * nu * self.m_uv) @ x
        return mins, mom, eqr, obj

    def x_from_zfulls(self, zf):
        """zf: dict family -> full hat tensor (23,)*n (site-symmetrised)"""
        A = self.A; offs = np.cumsum([0] + [f * f for f in A.f]); x = np.zeros(self.N); res = 0.0
        for c in FAMS:
            for s in self.fam[c]['srcs']:
                lam = s['lam']; n = len(lam)
                Zl = zf[c][np.ix_(*[np.arange(offs[l], offs[l + 1]) for l in lam])]
                Zl = Zl.reshape([A.f[l] for l in lam for _ in (0, 1)])
                Zl = Zl.transpose([2 * q for q in range(n)] + [2 * q + 1 for q in range(n)]).reshape(s['D'], s['D'])
                rec = np.zeros_like(Zl); k0 = 0
                for cl in s['classes']:
                    V = reduce_blocks(Zl[None], cl)[0]
                    rec += cl['U'] @ np.kron(V, np.eye(cl['d'])) @ cl['U'].T
                    vals = svec((V + V.T) / 2); x[s['off'] + k0: s['off'] + k0 + len(vals)] = vals; k0 += len(vals)
                res = max(res, np.abs(rec - Zl).max())
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
            # column index of replica r at site s is the row index of replica sigma_s^{-1}(r)
            subs.append(''.join(X[inv(sig[s])[r]][s] for s in range(n - 1)) + b[r] + R[r]); ops.append(Ut[r].conj())
        M = np.einsum(','.join(subs) + '->' + ''.join(a) + ''.join(b), *ops, optimize='greedy').reshape(81, 81)
        out[head] = PT_flat @ M.ravel()
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]); mu = float(sys.argv[2]); nu = float(sys.argv[3]); s = float(sys.argv[4])
    cons = sys.argv[5].split(',') if len(sys.argv) > 5 else ['SWAP', 'ORTH', 'MARG', 'PT', 'COV', 'LE2']
    algs = build_algebras()
    rel = Herm2(n, algs, cons)
    sol = rel.solve(mu, nu, s)
    x = np.array(sol.x)
    mins, mom, eqr, obj = rel.evaluate(x, mu, nu, s)
    print(f"herm2 n={n} C = {mu}|u><u| {'+' if s > 0 else '-'} {nu}|v><v|, cons {'+'.join(cons)}: status {sol.status}, iters {sol.iterations}, "
          f"q relaxation = {sol.obj_val:.10e} (build {rel.build_time:.1f}s, solve {rel.solve_time:.1f}s)")
    print("   min eig:", {k: f"{v:.1e}" for k, v in mins.items()}, f"moment res {mom:.1e}, eq res {eqr:.1e}")
