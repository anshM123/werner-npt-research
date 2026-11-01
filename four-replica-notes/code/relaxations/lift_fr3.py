"""Facially reduced lifted (all-ratio) certificate at n = 4 in (+, -) coordinates, floating point.

With w+ = (1, 1), w- = (1, -1) put Y+ = X11 + 2 X12 + X22, Y- = X11 - 2 X12 + X22, Z = X11 - X22 and
eps+ = w+^T E w+, eps- = w-^T E w-, zeta = w+^T E w-.  The lifted problem of lift2.py becomes
   E Y+ = e eps+,  E Y- = e eps-,  E Z = e zeta,  [[eps+, zeta], [zeta, eps-]] >= 0,  (eps+ + eps-)/2 = 1,
   block k:  [[N+^T W(Y+) N+, N+^T W(Z) N-], [N-^T W(Z) N+, N-^T W(Y-) N-]] >= 0      (after facial reduction)
   minimise  (1/4) [(m_uu + 2 m_uv + m_vv).Y+ + (m_uu - 2 m_uv + m_vv).Y- + 2 (m_uu - m_vv).Z]
with N+-_k the scalar faces at t = +-1.  Every block entry depends on one of Y+, Y-, Z only.  Each block is cut to
V+ (+) V-, V+ = sum of the ranges of N+^T W(e_i) N+ and N+^T W(e_i) N-, V- likewise (the data span, exactly
equivalent).  Solved with Clarabel (faer).  Zero points from lift_fr2_zeros.npz (made by lift_fr2.py).
usage: python lift_fr3.py [filter] [solver: faer|qdldl]"""
import sys, os, time, numpy as np, scipy.sparse as sp, scipy.linalg as sla, clarabel
os.environ['FR_REAL'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def rangebasis(mats_list, d, rng, tol=1e-10):
    """orthonormal basis of the sum of column spaces of the given d x * matrices (stacked along the last axis)"""
    acc = np.zeros((d, d + 4))
    for M in mats_list:                      # M: d x c x n  (n matrices of size d x c)
        if M.size == 0: continue
        for c0 in range(0, M.shape[2], 256):
            Mc = M[:, :, c0:c0 + 256]
            G = rng.standard_normal((Mc.shape[1], d + 4, Mc.shape[2]))
            acc += np.einsum('abi,bci->ac', Mc, G)
    if not np.any(acc): return np.zeros((d, 0))
    U, s, _ = np.linalg.svd(acc, full_matrices=False)
    return U[:, :int((s > tol * max(1.0, s[0])).sum())]

if __name__ == "__main__":
    filt = sys.argv[1] if len(sys.argv) > 1 else 'COV'
    method = sys.argv[2] if len(sys.argv) > 2 else 'faer'
    t0 = time.time()
    from alg4 import build_algebras
    from herm4 import Herm4, ALL_CONS
    from rep4 import unsvec, svec, svec_pairs
    import fr2_herm4 as FR
    keep = FR.FILTERS[filt]
    algs = build_algebras(); rel = Herm4(4, algs, ALL_CONS, verbose=False); N = rel.N
    dat = np.load('lift_fr2_zeros.npz'); xp = list(dat['xp']); xm = list(dat['xm'])
    print(f"n=4 lifted FR (+-) ({filt}, {method}): build {time.time() - t0:.0f}s; zeros {len(xp)} (t=+1), {len(xm)} (t=-1), "
          f"max |q| {max(abs(rel.evaluate(x, 1.0, 1.0)[3]) for x in xp):.1e} / {max(abs(rel.evaluate(x, 1.0, -1.0)[3]) for x in xm):.1e}",
          flush=True)
    xsp = sum(xp); xsm = sum(xm)
    NX = 3 * N + 3; EP, EM, ZE = 3 * N, 3 * N + 1, 3 * N + 2          # columns: Y+ | Y- | Z | eps+ eps- zeta
    Emat = sp.csr_matrix(np.vstack([rel.mom_rows, rel.Q])); evec = np.concatenate([rel.mom_rhs, np.zeros(rel.Q.shape[0])])
    ne = Emat.shape[0]
    rows = []; b = []; cones = []
    for off, ecol in ((0, EP), (N, EM), (2 * N, ZE)):
        R = sp.hstack([sp.csr_matrix((ne, off)), Emat, sp.csr_matrix((ne, NX - off - N))], format='lil')
        R[:, ecol] = -evec.reshape(-1, 1); rows.append(R.tocsr())
    rows.append(sp.csr_matrix((np.array([0.5, 0.5]), (np.array([0, 0]), np.array([EP, EM]))), shape=(1, NX)))
    b += [0.0] * (3 * ne) + [1.0]; cones.append(clarabel.ZeroConeT(3 * ne + 1))
    rows.append(sp.csr_matrix((np.array([-1.0, -np.sqrt(2), -1.0]), (np.arange(3), np.array([EP, ZE, EM]))), shape=(3, NX)))
    b += [0.0] * 3; cones.append(clarabel.PSDTriangleConeT(2))
    specs = [('X', o, m, None) for (o, m) in rel.wcones] + [(name, None, m, Md) for (name, mu, m, Md) in rel.blocks]
    sizes = []; info = []; tot0 = 0; rng = np.random.default_rng(11)
    for (name, o, m, Md) in specs:
        if name != 'X' and not keep(name): continue
        k = m * (m + 1) // 2; Mf = np.zeros((k, N))
        if name == 'X': Mf[np.arange(k), o + np.arange(k)] = 1.0
        else:
            for cf, Mc in Md.items(): Mf[:, rel.fam[cf]['start']:rel.fam[cf]['stop']] = Mc
        def ker(xs):
            ev, V = np.linalg.eigh(unsvec(Mf @ xs, m)); return V[:, ev < 1e-9 * max(1.0, ev.max())]
        Np, Nm = ker(xsp), ker(xsm); dp, dm = Np.shape[1], Nm.shape[1]
        if dp + dm == 0: continue
        nzc = np.nonzero(np.abs(Mf).max(axis=0) > 0)[0]               # parameters that enter this block
        pm = svec_pairs(m); I = np.array([i for i, j in pm]); J = np.array([j for i, j in pm])
        S = np.array([1.0 if i == j else 1.0 / np.sqrt(2.0) for i, j in pm])
        Wt = np.zeros((m, m, len(nzc))); Wt[I, J, :] = Mf[:, nzc] * S[:, None]; Wt[J, I, :] = Mf[:, nzc] * S[:, None]
        Dp = np.einsum('ap,abi,bq->pqi', Np, Wt, Np); Dm = np.einsum('ap,abi,bq->pqi', Nm, Wt, Nm)
        Gx = np.einsum('ap,abi,bq->pqi', Np, Wt, Nm)                   # d+ x d- x n
        Vp = rangebasis([Dp, Gx], dp, rng) if dp else np.zeros((0, 0))
        Vm = rangebasis([Dm, np.transpose(Gx, (1, 0, 2))], dm, rng) if dm else np.zeros((0, 0))
        rp, rm = Vp.shape[1], Vm.shape[1]; r = rp + rm; tot0 += dp + dm
        if r == 0: continue
        Dp = np.einsum('ap,abi,bq->pqi', Vp, Dp, Vp) if rp else None
        Dm = np.einsum('ap,abi,bq->pqi', Vm, Dm, Vm) if rm else None
        Gx = np.einsum('ap,abi,bq->pqi', Vp, Gx, Vm) if (rp and rm) else None
        pr = svec_pairs(r); R = np.zeros((len(pr), NX))
        for t_, (p, q) in enumerate(pr):
            if q < rp: R[t_, nzc] = -Dp[p, q, :] * (1.0 if p == q else np.sqrt(2.0))
            elif p >= rp: R[t_, N + nzc] = -Dm[p - rp, q - rp, :] * (1.0 if p == q else np.sqrt(2.0))
            else: R[t_, 2 * N + nzc] = -Gx[p, q - rp, :] * np.sqrt(2.0)
        R[np.abs(R) < 1e-13 * max(1.0, np.abs(R).max())] = 0.0
        rows.append(sp.csr_matrix(R)); b += [0.0] * R.shape[0]; cones.append(clarabel.PSDTriangleConeT(r))
        sizes.append(r); info.append((name, dp, dm, rp, rm))
        del Mf, Wt, R
    A = sp.vstack(rows, format='csc'); bvec = np.array(b); del rows
    print(f"   lifted reduced (+-): face dim {tot0} -> {sum(sizes)} after spans; {A.shape[0]} rows, nnz {A.nnz / 1e6:.1f}M, "
          f"{len(sizes)} PSD blocks (largest {max(sizes)}), sum k^2 {sum((s * (s + 1) // 2) ** 2 for s in sizes) / 1e6:.1f}M; "
          f"{time.time() - t0:.0f}s", flush=True)
    c = np.zeros(NX)
    c[0:N] = (rel.m_uu + 2 * rel.m_uv + rel.m_vv) / 4; c[N:2 * N] = (rel.m_uu - 2 * rel.m_uv + rel.m_vv) / 4
    c[2 * N:3 * N] = (rel.m_uu - rel.m_vv) / 2
    G = (A.T @ A).toarray(); scale = np.sqrt(np.maximum(np.diag(G), 1e-300)); Gs = G / np.outer(scale, scale)
    ev, V = np.linalg.eigh(Gs); null = ev < 1e-12 * ev[-1]; rank = int((~null).sum())
    cs = c / scale; cnull = np.linalg.norm(V[:, null].T @ cs) / np.linalg.norm(cs) if np.any(null) else 0.0
    Lr = V[:, ~null] * np.sqrt(ev[~null]); _, _, perm = sla.qr(Lr.T, pivoting=True, mode='economic'); Jc = np.sort(perm[:rank])
    Am = A[:, Jc]; cobj = c[Jc]; del G, Gs, V, Lr
    print(f"   rank {rank} of {NX} columns; relative |c on null| {cnull:.1e}", flush=True)
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 500
    st.tol_gap_abs = 1e-9; st.tol_gap_rel = 1e-9; st.tol_feas = 1e-9
    st.chordal_decomposition_enable = False; st.direct_solve_method = method
    t1 = time.time()
    sol = clarabel.DefaultSolver(sp.csc_matrix((len(Jc), len(Jc))), cobj, Am, bvec, cones, st).solve()
    z = np.array(sol.z); res = Am.T @ z + cobj
    row = 3 * ne + 1; zE = z[row:row + 3]; row += 3; mz = np.inf; worst = None; small = 0; eigs = []
    mpp = np.inf; mmm = np.inf
    for s, nm in zip(sizes, info):
        kk = s * (s + 1) // 2; Sb = unsvec(z[row:row + kk], s); row += kk
        e = np.linalg.eigvalsh(Sb).min(); eigs.append(e)
        rp = nm[3]
        if rp: mpp = min(mpp, np.linalg.eigvalsh(Sb[:rp, :rp]).min())
        if s - rp: mmm = min(mmm, np.linalg.eigvalsh(Sb[rp:, rp:]).min())
        if e < 1e-6: small += 1
        if e < mz: mz, worst = e, nm
    print(f"   lifted FR (+-) ({filt}): status {sol.status}, iters {sol.iterations}, primal {sol.obj_val:+.3e}, "
          f"dual {-(bvec @ z):+.3e}, ||A^T z + c|| {np.abs(res).max():.1e}; min eig lifted dual blocks {mz:.2e} (at {worst}); "
          f"blocks with eig < 1e-6: {small} of {len(sizes)}; min eig ++ parts {mpp:.2e}, -- parts {mmm:.2e}; dual E-block "
          f"{np.round(zE, 5).tolist()}; solve {time.time() - t1:.0f}s", flush=True)
    np.save('lift_fr3_eigs.npy', np.array(eigs))
