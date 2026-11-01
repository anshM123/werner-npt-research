"""Certificate for q_4(a P_u + b P_v) >= 0 for all real a, b (lifted problem in (+, -) coordinates) from the exact
lifted, face- and span-reduced data (facedataLr.pkl) and the exact linear data (linL.pkl).

Certificate: y+, y-, yZ (rational, one entry per equality) and lifted dual blocks S_k = [[A_k, B_k], [B_k^T, C_k]] >= 0 with
    c+/4 = E^T y+ + sum_k F+_k^*(A_k),   c-/4 = E^T y- + sum_k F-_k^*(C_k),   d/2 = E^T yZ + sum_k 2 G_k^*(B_k),
    e . y+ = e . y- = e . yZ = 0,
F+_k^*(A)_i = Tr(A F+_k(e_i)), G_k^*(B)_i = Tr(B^T G_k(e_i)).  Then on the lifted relaxation
    (1/4)(c+ . Y+ + c- . Y- + 2 d . Z) = (e . y+) eps+ + (e . y-) eps- + (e . yZ) zeta + sum_k Tr(S_k L_k) = sum_k Tr(S_k L_k) >= 0.
Steps: (1) float lifted problem (variables Y+, Y-, Z, eps+, eps-, zeta; equalities, the 2 x 2 eps-block, the
normalisation (eps+ + eps-)/2 = 1, all lifted blocks with Jacobi scaling), independent columns by pivoted QR (null
directions dropped), Clarabel with the faer direct solver, chordal decomposition off; (2) float polish of the dual:
minimum-norm correction of the dual residual over ALL columns (pseudo-inverse of the scaled Gram matrix); (3) rounding
of y and S_k to dyadic grids; (4) exact correction, separately for the three decoupled identities (+, -, Z), on column
subsets of full row rank (pivoted QR, exact square solve with flint Dixon); (5) exact check of all 3 x 1332 equations;
(6) exact LDL^T of every S_k; (7) check at an exact entangled true point for a few ratios.  Saves certL.pkl.
usage: python certL.py [--tol 1e-10] [--mode whiten|colscale|plain]   (default whiten)"""
import sys, os, time, pickle, math, numpy as np, flint, scipy.sparse as sp, scipy.linalg, clarabel
from gmpy2 import mpq

def log(s): print(s, flush=True)

def svec_pairs(m): return [(i, j) for j in range(m) for i in range(j + 1)]
def unsvec(v, m):
    V = np.zeros((m, m)); k = 0
    for j in range(m):
        for i in range(j + 1):
            if i == j: V[i, i] = v[k]
            else: V[i, j] = V[j, i] = v[k] / np.sqrt(2.0)
            k += 1
    return V

def exact_ldl(S):
    """exact LDL^T of a symmetric rational matrix (list of lists of mpq); returns (psd, smallest positive pivot, #zero pivots)"""
    k = len(S); M = [row[:] for row in S]; minpiv = None; nzero = 0
    for i in range(k):
        p = M[i][i]
        if p < 0: return False, p, nzero
        if p == 0:
            if any(M[i][j] != 0 for j in range(i + 1, k)): return False, p, nzero
            nzero += 1; continue
        minpiv = p if minpiv is None or p < minpiv else minpiv
        for j in range(i + 1, k):
            if M[j][i] != 0:
                f = M[j][i] / p
                for l in range(i + 1, k): M[j][l] -= f * M[i][l]
    return True, minpiv, nzero

def main():
    t0 = time.time()
    tol = 1e-10; mode = 'whiten'
    if '--tol' in sys.argv: tol = float(sys.argv[sys.argv.index('--tol') + 1])
    if '--mode' in sys.argv: mode = sys.argv[sys.argv.index('--mode') + 1]
    assert mode in ('plain', 'colscale', 'whiten')
    lin = pickle.load(open('linL.pkl', 'rb')); fd = pickle.load(open('facedataLr.pkl', 'rb'))
    N = lin['N']; E = lin['E']; e = lin['e']; neq = len(E)
    tgt = {'+': [v / 4 for v in lin['cplus']], '-': [v / 4 for v in lin['cminus']], 'Z': [v / 2 for v in lin['d']]}
    blocks = fd['blocks']; nb = len(blocks)
    sizes = [b['rp'] + b['rm'] for b in blocks]
    log(f"[0] data: {N} parameters, {neq} equalities, {nb} lifted blocks ({sum(1 for b in blocks if b['kind'] == 'COV')} COV, "
        f"{sum(1 for b in blocks if b['kind'] == 'X')} X), total dimension {sum(sizes)}, largest {max(sizes)}; {time.time() - t0:.0f}s")
    # ------------------------------------------------------------ float data
    bf = []
    for b in blocks:
        out = []
        for gi, idx, Fp, Fm, G, D in b['data']:
            f = lambda M: (np.array(M, dtype=float) / float(D)) if M.size else np.zeros(M.shape)
            out.append((np.array(idx), f(Fp), f(Fm), f(G)))
        bf.append(out)
    Ef = np.array([[float(v) for v in row] for row in E]); ef = np.array([float(v) for v in e])
    esc = np.linalg.norm(Ef, axis=1); Efs = Ef / esc[:, None]; efs = ef / esc
    NX = 3 * N + 3; EP, EM, ZE = 3 * N, 3 * N + 1, 3 * N + 2
    OFF = {'+': 0, '-': N, 'Z': 2 * N}
    cvec = np.zeros(NX)
    for g in ('+', '-', 'Z'): cvec[OFF[g]:OFF[g] + N] = [float(v) for v in tgt[g]]
    # zero cone: E Y+ - e eps+ = 0, E Y- - e eps- = 0, E Z - e zeta = 0, (eps+ + eps-)/2 = 1
    R0 = []; C0 = []; V0 = []
    for blk, (off, ecol) in enumerate(((0, EP), (N, EM), (2 * N, ZE))):
        jj, ii = np.nonzero(Efs)
        R0 += list(blk * neq + jj); C0 += list(off + ii); V0 += list(Efs[jj, ii])
        R0 += list(blk * neq + np.arange(neq)); C0 += [ecol] * neq; V0 += list(-efs)
    R0 += [3 * neq, 3 * neq]; C0 += [EP, EM]; V0 += [0.5, 0.5]
    nz0 = 3 * neq + 1
    bvec = [0.0] * (3 * neq) + [1.0]
    # eps block (PSD 2 x 2): svec([[eps+, zeta], [zeta, eps-]]) = (eps+, sqrt2 zeta, eps-)
    R0 += [nz0, nz0 + 1, nz0 + 2]; C0 += [EP, ZE, EM]; V0 += [-1.0, -np.sqrt(2.0), -1.0]
    bvec += [0.0] * 3
    rowc = nz0 + 3
    # lifted blocks: scalar ones first (one nonnegative cone), then PSD ones
    order = [k for k in range(nb) if sizes[k] == 1] + [k for k in range(nb) if sizes[k] > 1]
    trans = {}; Rb = []; Cb = []; Vb = []; rowstart = {}
    for k in order:
        b = blocks[k]; rp, rm = b['rp'], b['rm']; r = rp + rm
        dg2 = np.zeros(r)
        for idx, Fp, Fm, G in bf[k]:
            if rp: dg2[:rp] += np.einsum('kab,kab->a', Fp, Fp) + (np.einsum('kab,kab->a', G, G) if rm else 0)
            if rm: dg2[rp:] += np.einsum('kab,kab->a', Fm, Fm) + (np.einsum('kab,kab->b', G, G) if rp else 0)
        assert np.all(dg2 > 0), f"block {b['name']} has a vacuous direction"
        dg = np.sqrt(dg2)
        P = svec_pairs(r); I = np.array([p[0] for p in P]); J = np.array([p[1] for p in P])
        fac = np.where(I == J, 1.0, np.sqrt(2.0)) / (dg[I] * dg[J])
        pp = np.nonzero(J < rp)[0]; mm = np.nonzero(I >= rp)[0]; pm = np.nonzero((I < rp) & (J >= rp))[0]
        rr = []; cc = []; vv = []
        for idx, Fp, Fm, G in bf[k]:
            K = len(idx)
            for sel, M, off, ia, ja in ((pp, Fp, 0, I, J), (mm, Fm, N, I - rp, J - rp), (pm, G, 2 * N, I, J - rp)):
                if len(sel) == 0 or M.size == 0: continue
                vals = M[:, ia[sel], ja[sel]] * fac[sel][None, :]          # K x len(sel)
                kk, ss = np.nonzero(vals)
                rr.append(sel[ss]); cc.append(off + idx[kk]); vv.append(vals[kk, ss])
        rr = np.concatenate(rr); cc = np.concatenate(cc); vv = np.concatenate(vv)
        s_k = np.abs(vv).max(); vv = vv / s_k
        trans[k] = (dg, s_k); rowstart[k] = rowc
        Rb.append(rowc + rr); Cb.append(cc); Vb.append(-vv)
        rowc += len(P); bvec += [0.0] * len(P)
    A = sp.csc_matrix((np.concatenate([np.array(V0)] + Vb), (np.concatenate([np.array(R0)] + Rb), np.concatenate([np.array(C0)] + Cb))),
                      shape=(rowc, NX))
    bvec = np.array(bvec)
    nsc = sum(1 for k in order if sizes[k] == 1)
    cones = [clarabel.ZeroConeT(nz0), clarabel.PSDTriangleConeT(2)]
    if nsc: cones.append(clarabel.NonnegativeConeT(nsc))
    cones += [clarabel.PSDTriangleConeT(sizes[k]) for k in order if sizes[k] > 1]
    log(f"[1] float lifted problem: {A.shape[0]} rows x {NX} columns, nnz {A.nnz / 1e6:.1f}M; {nsc} scalar blocks, "
        f"{nb - nsc} PSD blocks (largest {max(sizes)}), sum svec^2 {sum((s * (s + 1) // 2) ** 2 for s in sizes) / 1e6:.1f}M; "
        f"{time.time() - t0:.0f}s")
    # ------------------------------------------------------------ independent columns
    # Gram matrix A^T A accumulated block by block with dense products (each block touches few columns)
    Gm = np.zeros((NX, NX))
    A0 = sp.csc_matrix((np.array(V0), (np.array(R0), np.array(C0))), shape=(nz0 + 3, NX)).toarray()
    Gm += A0.T @ A0; del A0
    for rr_, cc_, vv_ in zip(Rb, Cb, Vb):
        ucols, loc = np.unique(cc_, return_inverse=True)
        r0 = rr_.min(); Ak = np.zeros((rr_.max() - r0 + 1, len(ucols))); np.add.at(Ak, (rr_ - r0, loc), vv_)
        Gm[np.ix_(ucols, ucols)] += Ak.T @ Ak
    del Rb, Cb, Vb
    scale = np.sqrt(np.maximum(np.diag(Gm), 1e-300)); Gs = Gm / np.outer(scale, scale); del Gm
    log(f"    Gram matrix built; {time.time() - t0:.0f}s")
    ev, V = np.linalg.eigh(Gs); null = ev < 1e-12 * ev[-1]; rank = int((~null).sum())
    cs = cvec / scale
    cnull = np.linalg.norm(V[:, null].T @ cs) / np.linalg.norm(cs) if np.any(null) else 0.0
    ngrp = {g: int(np.sum(np.abs(V[OFF[g]:OFF[g] + N, null]).max(axis=0) > 1e-8)) for g in OFF} if np.any(null) else {}
    Lr = V[:, ~null] * np.sqrt(ev[~null]); _, _, perm = scipy.linalg.qr(Lr.T, pivoting=True, mode='economic')
    Jc = np.sort(perm[:rank]); del Lr
    log(f"    rank {rank} of {NX} columns ({NX - rank} null directions; null vectors touching Y+ / Y- / Z: "
        f"{ngrp.get('+', 0)} / {ngrp.get('-', 0)} / {ngrp.get('Z', 0)}); relative objective component on the null space "
        f"{cnull:.1e}; column norms {scale.min():.1e} .. {scale.max():.1e}; |c_j| / norm_j up to {np.abs(cs).max():.1e}; "
        f"scaled Gram eigenvalues (non-null) {ev[~null].min() / ev[-1]:.1e} .. 1; {time.time() - t0:.0f}s")
    st = clarabel.DefaultSettings(); st.verbose = True; st.max_iter = 500
    st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
    st.chordal_decomposition_enable = False; st.direct_solve_method = 'faer'
    t1 = time.time()
    if mode == 'plain':                       # selected columns as they are
        Amat = A[:, Jc].tocsc(); qv = cvec[Jc]
    elif mode == 'colscale':                  # selected columns scaled to unit norm (same dual z)
        Amat = (A[:, Jc] @ sp.diags(1.0 / scale[Jc])).tocsc(); qv = cvec[Jc] / scale[Jc]
    else:                                     # 'whiten': per column group (Y+, Y-, Z) whitening, eps columns scaled
        parts = []; qs = []
        for g in ('+', '-', 'Z'):
            cols = OFF[g] + np.arange(N)
            evg, Vg = np.linalg.eigh(Gs[np.ix_(cols, cols)]); keep = evg > 1e-12 * evg[-1]
            Wg = (Vg[:, keep] / np.sqrt(evg[keep])) / scale[cols][:, None]
            Ag = A[:, cols]; rows_g = np.unique(Ag.nonzero()[0])
            Dg = np.zeros((A.shape[0], Wg.shape[1])); Dg[rows_g] = Ag[rows_g] @ Wg
            Dg[np.abs(Dg) < 1e-15 * np.abs(Dg).max()] = 0.0
            parts.append(sp.csc_matrix(Dg)); qs.append(Wg.T @ cvec[cols]); del Dg
            log(f"    whitening {g}: rank {int(keep.sum())} of {N}")
        ecols = np.array([EP, EM, ZE])
        parts.append((A[:, ecols] @ sp.diags(1.0 / scale[ecols])).tocsc()); qs.append(cvec[ecols] / scale[ecols])
        Amat = sp.hstack(parts, format='csc'); qv = np.concatenate(qs); del parts
    log(f"    solver input ({mode}): {Amat.shape[0]} x {Amat.shape[1]}, nnz {Amat.nnz / 1e6:.1f}M; {time.time() - t0:.0f}s")
    sol = clarabel.DefaultSolver(sp.csc_matrix((Amat.shape[1], Amat.shape[1])), qv, Amat, bvec, cones, st).solve()
    del Amat
    z = np.array(sol.z)
    res_sel = np.abs(A[:, Jc].T @ z + cvec[Jc]).max(); res_all = np.abs(A.T @ z + cvec).max()
    log(f"[2] float solve: status {sol.status}, iterations {sol.iterations}, primal {sol.obj_val:+.3e}, dual {-(bvec @ z):+.3e}, "
        f"||A^T z + c|| {res_sel:.1e} (selected columns), {res_all:.1e} (all columns); solve {time.time() - t1:.0f}s")
    # ------------------------------------------------------------ float polish (minimum-norm, all columns, scaled)
    def dual_stats(z):
        mins = []; minpp = np.inf; minmm = np.inf; worst = None
        for k in order:
            r = sizes[k]; kk = r * (r + 1) // 2; Sb = unsvec(z[rowstart[k]:rowstart[k] + kk], r)
            m_ = np.linalg.eigvalsh(Sb).min(); mins.append(m_)
            rp = blocks[k]['rp']
            if rp: minpp = min(minpp, np.linalg.eigvalsh(Sb[:rp, :rp]).min())
            if r - rp: minmm = min(minmm, np.linalg.eigvalsh(Sb[rp:, rp:]).min())
            if worst is None or m_ < worst[0]: worst = (m_, blocks[k]['name'], blocks[k]['rp'], blocks[k]['rm'])
        return min(mins), worst, minpp, minmm, sum(1 for m_ in mins if m_ < 1e-6)
    mz, worst, mpp, mmm, nsmall = dual_stats(z)
    zE = z[nz0:nz0 + 3]
    log(f"    smallest lifted dual eigenvalue (Jacobi-normalised blocks) {mz:.3e} at {worst[1]} (r+ {worst[2]}, r- {worst[3]}); "
        f"blocks with eigenvalue < 1e-6: {nsmall}; smallest ++ part {mpp:.3e}, -- part {mmm:.3e}; eps-block dual "
        f"{np.round(zE, 8).tolist()}, normalisation multiplier {z[3 * neq]:+.2e}")
    for it in range(2):
        rr_ = -(A.T @ z + cvec); rs = rr_ / scale
        w = V[:, ~null].T @ rs; w = w / ev[~null]; dz = (A @ ((V[:, ~null] @ w) / scale))
        z = z + dz
        log(f"    polish {it + 1}: max |dz| {np.abs(dz).max():.1e}; ||A^T z + c|| now {np.abs(A.T @ z + cvec).max():.1e}")
    mz2, worst2, _, _, _ = dual_stats(z)
    log(f"    after polish: smallest lifted dual eigenvalue {mz2:.3e}; {time.time() - t0:.0f}s")
    if mz2 <= 1e-9:
        log("n=4 lifted: FAILED (lifted dual not strictly feasible: face not minimal or solve inaccurate)"); return
    yf = {'+': -z[0:neq] / esc, '-': -z[neq:2 * neq] / esc, 'Z': -z[2 * neq:3 * neq] / esc}
    Sf = {}
    for k in order:
        r = sizes[k]; kk = r * (r + 1) // 2; dg, s_k = trans[k]
        Sb = unsvec(z[rowstart[k]:rowstart[k] + kk], r)
        Sf[k] = Sb / np.outer(dg, dg) / s_k
    with open('float_dualL.pkl', 'wb') as fh: pickle.dump(dict(y=yf, S=Sf, mineig=mz2), fh)
    # ------------------------------------------------------------ rounding to dyadic grids
    def rq(v, ex):
        if ex >= 0: return mpq(int(round(float(v) * 2.0 ** ex)), 2 ** ex)
        return mpq(int(round(float(v) / 2.0 ** (-ex))) * 2 ** (-ex))
    yq = {}
    for g in yf:
        assert np.all(np.isfinite(yf[g]))
        ex = 40 - int(np.floor(np.log2(np.abs(yf[g]).max() + 1e-300))); yq[g] = [rq(v, ex) for v in yf[g]]
    Sq = {}
    for k, S in Sf.items():
        assert np.all(np.isfinite(S))
        ex = 40 - int(np.floor(np.log2(np.abs(S).max() + 1e-300)))
        r = S.shape[0]; M = [[None] * r for _ in range(r)]
        for a_ in range(r):
            for b_ in range(a_, r): M[a_][b_] = M[b_][a_] = rq(S[a_, b_], ex)
        Sq[k] = M
    # ------------------------------------------------------------ exact adjoints and residuals
    def part(k, g, M):
        """integer numerator array of the relevant part of S_k and its common denominator"""
        rp = blocks[k]['rp']; r = sizes[k]
        if g == '+': rows, cols = range(0, rp), range(0, rp)
        elif g == '-': rows, cols = range(rp, r), range(rp, r)
        else: rows, cols = range(0, rp), range(rp, r)
        den = 1
        for a_ in rows:
            for b_ in cols: den = den * int(M[a_][b_].denominator) // math.gcd(den, int(M[a_][b_].denominator))
        arr = np.array([[int(M[a_][b_] * den) for b_ in cols] for a_ in rows], dtype=object).reshape(len(rows), len(cols))
        return arr, den
    def adjoint(k, g, M):
        """F+^*(A) (g '+'), F-^*(C) (g '-'), 2 G^*(B) (g 'Z') as dict i -> mpq"""
        out = {}; arr, den = part(k, g, M)
        if arr.size == 0: return out
        for gi, idx, Fp, Fm, G, D in blocks[k]['data']:
            Md = {'+': Fp, '-': Fm, 'Z': G}[g]
            if Md.size == 0: continue
            vals = Md.reshape(Md.shape[0], -1).dot(arr.ravel())
            f = 2 if g == 'Z' else 1
            for kk, i in enumerate(idx):
                if vals[kk] != 0: out[i] = out.get(i, mpq(0)) + mpq(int(vals[kk]) * f, D * den)
        return out
    def residual(g):
        rres = [mpq(v) for v in tgt[g]]
        for j in range(neq):
            yj = yq[g][j]
            if yj == 0: continue
            for i, v in enumerate(E[j]):
                if v != 0: rres[i] -= v * yj
        for k in Sq:
            for i, v in adjoint(k, g, Sq[k]).items(): rres[i] -= v
        robj = -sum((e[j] * yq[g][j] for j in range(neq)), mpq(0))
        return rres, robj
    for g in ('+', '-', 'Z'):
        rres, robj = residual(g)
        log(f"[3] rounded dual, identity {g}: residual {max(abs(float(v)) for v in rres):.1e}, e . y residual {float(robj):.1e}")
    # ------------------------------------------------------------ exact correction, per identity
    corr_info = {}
    for g in ('+', '-', 'Z'):
        rres, robj = residual(g)
        cols = []; Fc = []
        for j in range(neq):
            v = np.concatenate([Ef[j], [ef[j]]]); cols.append(('y', j)); Fc.append(v / np.linalg.norm(v))
        for k in order:
            rp = blocks[k]['rp']; rm = blocks[k]['rm']
            if g == '+': prs = [(a_, b_) for (a_, b_) in svec_pairs(rp)]
            elif g == '-': prs = [(rp + a_, rp + b_) for (a_, b_) in svec_pairs(rm)]
            else: prs = [(a_, rp + b_) for a_ in range(rp) for b_ in range(rm)]
            if not prs: continue
            vk = np.zeros((len(prs), N + 1))
            for idx, Fp, Fm, G in bf[k]:
                for t_, (a_, b_) in enumerate(prs):
                    if g == '+': vk[t_, idx] += Fp[:, a_, b_] * (1.0 if a_ == b_ else 2.0)
                    elif g == '-': vk[t_, idx] += Fm[:, a_ - rp, b_ - rp] * (1.0 if a_ == b_ else 2.0)
                    else: vk[t_, idx] += 2.0 * G[:, a_, b_ - rp]
            for t_, (a_, b_) in enumerate(prs):
                nv = np.linalg.norm(vk[t_])
                if nv == 0: continue
                cols.append(('S', k, a_, b_)); Fc.append(vk[t_] / nv)
        Fm_ = np.array(Fc).T; del Fc
        Qf, Rf, piv = scipy.linalg.qr(Fm_, pivoting=True, mode='economic')
        dR = np.abs(np.diag(Rf)); rank_g = int(np.sum(dR > 1e-10 * dR[0]))
        chosen = [cols[j] for j in piv[:rank_g]]
        del Fm_, Qf, Rf
        def exact_col(cd):
            v = [mpq(0)] * (N + 1)
            if cd[0] == 'y':
                for i in range(N): v[i] = E[cd[1]][i]
                v[N] = e[cd[1]]
            else:
                _, k, a_, b_ = cd; rp = blocks[k]['rp']
                for gi, idx, Fp, Fm, G, D in blocks[k]['data']:
                    for kk, i in enumerate(idx):
                        if g == '+': val = Fp[kk, a_, b_] * (1 if a_ == b_ else 2)
                        elif g == '-': val = Fm[kk, a_ - rp, b_ - rp] * (1 if a_ == b_ else 2)
                        else: val = 2 * G[kk, a_, b_ - rp]
                        if val: v[i] += mpq(int(val), D)
            return v
        C = [exact_col(cd) for cd in chosen]
        rvec = rres + [robj]
        Cf = np.array([[float(x) for x in col] for col in C])
        _, _, rpiv = scipy.linalg.qr(Cf, pivoting=True, mode='economic')
        rowsel = sorted(int(j) for j in rpiv[:rank_g])
        cden = []
        for col in C:
            dd = 1
            for x in col: dd = dd * int(x.denominator) // math.gcd(dd, int(x.denominator))
            cden.append(dd)
        Asq = flint.fmpz_mat(rank_g, rank_g, [int(C[jj][ii] * cden[jj]) for ii in rowsel for jj in range(rank_g)])
        bsq = flint.fmpq_mat(rank_g, 1, [flint.fmpq(int(rvec[ii].numerator), int(rvec[ii].denominator)) for ii in rowsel])
        t1 = time.time()
        sol_d = flint.fmpq_mat(Asq).solve(bsq, algorithm='dixon')
        dl = [mpq(int(sol_d[jj, 0].p), int(sol_d[jj, 0].q)) * cden[jj] for jj in range(rank_g)]
        for cd, dv in zip(chosen, dl):
            if cd[0] == 'y': yq[g][cd[1]] += dv
            else:
                _, k, a_, b_ = cd; Sq[k][a_][b_] += dv
                if a_ != b_: Sq[k][b_][a_] += dv
        cbits = max(max(int(abs(v.numerator)).bit_length(), int(v.denominator).bit_length()) for v in dl)
        corr_info[g] = (rank_g, len(cols), max(abs(float(v)) for v in dl), cbits)
        log(f"[4] identity {g}: {len(cols)} columns, equation rank {rank_g} of {N + 1} (diag R ratio at cut "
            f"{dR[rank_g - 1] / dR[0]:.1e}, next {(dR[rank_g] / dR[0]) if rank_g < len(dR) else 0:.1e}); "
            f"{sum(1 for cd in chosen if cd[0] == 'y')} y's, {sum(1 for cd in chosen if cd[0] == 'S')} dual entries; exact square "
            f"solve {rank_g} x {rank_g} in {time.time() - t1:.0f}s; max |delta| {corr_info[g][2]:.1e}, correction bits up to {cbits}; "
            f"{time.time() - t0:.0f}s")
    exact_ok = True
    for g in ('+', '-', 'Z'):
        rres, robj = residual(g)
        ok = all(v == 0 for v in rres) and robj == 0
        exact_ok = exact_ok and ok
        log(f"[5] identity {g} exact in all {N} coordinates and e . y{'' if g != 'Z' else 'Z'} = 0: {ok}")
    with open('certL.pkl', 'wb') as fh:
        pickle.dump(dict(y={g: [(int(v.numerator), int(v.denominator)) for v in yq[g]] for g in yq},
                         S={k: [[(int(v.numerator), int(v.denominator)) for v in row] for row in M] for k, M in Sq.items()}), fh)
    if not exact_ok:
        log("n=4 lifted: FAILED (identities not exact after the correction)"); return
    # ------------------------------------------------------------ exact LDL^T
    allpsd = True; minpiv = None; t1 = time.time(); nzp = 0
    for k, M in Sq.items():
        ok, mp, nz = exact_ldl(M); nzp += nz
        if not ok: allpsd = False; log(f"    block {blocks[k]['name']}: NOT PSD"); break
        if mp is not None and (minpiv is None or mp < minpiv): minpiv = mp
    nbits = max(max(int(abs(v.numerator)).bit_length() for g in yq for v in yq[g]),
                max(int(abs(x.numerator)).bit_length() for M in Sq.values() for rw in M for x in rw))
    dbits = max(max(int(v.denominator).bit_length() for g in yq for v in yq[g]),
                max(int(x.denominator).bit_length() for M in Sq.values() for rw in M for x in rw))
    log(f"[6] exact LDL^T: all {len(Sq)} lifted dual blocks PSD {allpsd} (smallest pivot {float(minpiv) if minpiv else 0:.2e}, "
        f"zero pivots {nzp}); rational sizes up to {nbits} / {dbits} bits; LDL {time.time() - t1:.0f}s; {time.time() - t0:.0f}s")
    # ------------------------------------------------------------ check at an exact entangled true point
    rec = pickle.load(open('../exact4/pts/p_ent.pkl', 'rb'))
    xP = [mpq(p_, q_) for (p_, q_) in rec['x']]
    Ex = [sum((E[j][i] * xP[i] for i in range(N) if E[j][i] != 0), mpq(0)) for j in range(neq)]
    eqok = all(Ex[j] == e[j] for j in range(neq))
    dot = lambda c: sum((c[i] * xP[i] for i in range(N) if c[i] != 0), mpq(0))
    cpx, cmx, dx = dot(lin['cplus']), dot(lin['cminus']), dot(lin['d'])
    allok = eqok and cpx == mpq(*rec['q4'])
    adj = {g: {k: adjoint(k, g, Sq[k]) for k in Sq} for g in ('+', '-', 'Z')}
    tv = {g: {k: sum((v * xP[i] for i, v in adj[g][k].items()), mpq(0)) for k in Sq} for g in adj}
    for (a, b) in ((1, 0), (1, 1), (1, -1), (2, 1)):
        sp_, sm_, sz_ = (a + b) ** 2, (a - b) ** 2, a * a - b * b
        obj = (sp_ * cpx + sm_ * cmx + 2 * sz_ * dx) / 4
        terms = [sp_ * tv['+'][k] + sm_ * tv['-'][k] + sz_ * tv['Z'][k] for k in Sq]
        yE = sum(((sp_ * yq['+'][j] + sm_ * yq['-'][j] + sz_ * yq['Z'][j]) * Ex[j] for j in range(neq)), mpq(0))
        ok = obj == yE + sum(terms) and yE == 0 and all(t_ >= 0 for t_ in terms)
        allok = allok and ok
        log(f"[7] p_ent, ratio ({a},{b}): lifted objective {float(obj):.10f}; identity at the lifted point "
            f"{obj == yE + sum(terms)}; y-terms {yE}; all {len(terms)} terms >= 0: {all(t_ >= 0 for t_ in terms)} "
            f"(smallest {float(min(terms)):.2e})")
    if allpsd and exact_ok and allok:
        log(f"n=4 lifted: exact identities (+, -, Z) with e.y = 0, {len(Sq)} PSD lifted blocks (exact LDL^T), max {nbits}/{dbits} bits, "
            f"entangled true-point check passed -> CERTIFIED")
    else:
        log("n=4 lifted: FAILED (" + ("block not PSD" if not allpsd else "true-point check failed") + ")")

if __name__ == "__main__":
    main()
