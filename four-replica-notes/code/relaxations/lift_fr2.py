"""Facially reduced lifted (all-ratio) certificate at n = 4, floating point, memory-lean version of lift_fr.py.

Same lifted face as lift_fr.py: block k is compressed by B_k = [w+ (x) N+_k, w- (x) N-_k], N+-_k the scalar faces at
t = +-1 (kernels at the summed zeros).  In addition each compressed block is cut to the span V_k of its data
matrices (the union of the ranges of all D_i = data of parameter i): every primal block is supported on V_k x V_k, so
the dual block only matters on V_k (as in exact4/span4.py).  Zero points are cached in lift_fr2_zeros.npz.
usage: python lift_fr2.py [filter] [n_spinflip] [workers]"""
import sys, os, time, numpy as np, scipy.sparse as sp, scipy.linalg as sla, clarabel
os.environ['FR_REAL'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
from lift_fr import _zf_worker, pairs_tp1, pairs_tm1, svec_pair_map

def data_span(R, d, unsvec_idx, tol=1e-10, rng=None):
    """orthonormal basis of the sum of ranges of the symmetric d x d matrices whose svec are the columns of R"""
    cols = np.nonzero(np.abs(R).max(axis=0) > 0)[0]
    if len(cols) == 0: return np.zeros((d, 0))
    I, J, S = unsvec_idx
    acc = np.zeros((d, d + 4))
    for c0 in range(0, len(cols), 256):
        cc = cols[c0:c0 + 256]
        D = np.zeros((d, d, len(cc)))
        D[I, J, :] = R[:, cc] * S[:, None]; D[J, I, :] = R[:, cc] * S[:, None]
        G = rng.standard_normal((d, d + 4, len(cc)))
        acc += np.einsum('abi,bci->ac', D, G)
    U, s, _ = np.linalg.svd(acc, full_matrices=False)
    r = int((s > tol * max(1.0, s[0])).sum())
    return U[:, :r]

def unsvec_index(d, svec_pairs):
    pm = svec_pairs(d)
    I = np.array([i for (i, j) in pm]); J = np.array([j for (i, j) in pm])
    S = np.array([1.0 if i == j else 1.0 / np.sqrt(2.0) for (i, j) in pm])
    return I, J, S

if __name__ == "__main__":
    filt = sys.argv[1] if len(sys.argv) > 1 else 'COV'
    nsf = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    t0 = time.time(); cache = 'lift_fr2_zeros.npz'
    rng = np.random.default_rng(4242)
    PP = pairs_tp1(rng); PM = pairs_tm1(rng, nsf)
    job = None
    if not os.path.exists(cache):
        pool = Pool(nw); job = pool.map_async(_zf_worker, PP + PM)
    from alg4 import build_algebras
    from herm4 import Herm4, ALL_CONS
    from validate import F1_matrix
    from rep4 import unsvec, svec, svec_pairs
    import fr2_herm4 as FR
    keep = FR.FILTERS[filt]
    algs = build_algebras(); F1 = F1_matrix(algs[0])
    rel = Herm4(4, algs, ALL_CONS, verbose=False); N = rel.N
    print(f"n=4 lifted FR, span-reduced ({filt}): build {time.time() - t0:.0f}s ({N} parameters)", flush=True)
    if job is not None:
        xp = []; xm = []
        for k in (1, 3):
            for j in range(2): xp.append(rel.x_from_zfulls(FR.zero_point_k(4, k, rng, F1))[0])
        for k in (2, 4):
            for j in range(2): xm.append(rel.x_from_zfulls(FR.zero_point_k(4, k, rng, F1))[0])
        zfs = job.get(); pool.close()
        for i, zf in enumerate(zfs):
            (xp if i < len(PP) else xm).append(rel.x_from_zfulls(zf)[0])
        np.savez(cache, xp=np.array(xp), xm=np.array(xm))
    else:
        dat = np.load(cache); xp = list(dat['xp']); xm = list(dat['xm'])
    print(f"   zeros: {len(xp)} at t=+1 (max |q| {max(abs(rel.evaluate(x, 1.0, 1.0)[3]) for x in xp):.1e}), {len(xm)} at "
          f"t=-1 (max |q| {max(abs(rel.evaluate(x, 1.0, -1.0)[3]) for x in xm):.1e}); {time.time() - t0:.0f}s", flush=True)
    xsp = sum(xp); xsm = sum(xm)
    wp = np.array([1.0, 1.0]) / np.sqrt(2); wm = np.array([1.0, -1.0]) / np.sqrt(2)
    NX = 3 * N + 3; E11, E12, E22 = 3 * N, 3 * N + 1, 3 * N + 2
    Emat = sp.csr_matrix(np.vstack([rel.mom_rows, rel.Q])); evec = np.concatenate([rel.mom_rhs, np.zeros(rel.Q.shape[0])])
    ne = Emat.shape[0]
    blocks = []; b = []; cones = []
    for off, ecol in ((0, E11), (N, E12), (2 * N, E22)):
        R = sp.hstack([sp.csr_matrix((ne, off)), Emat, sp.csr_matrix((ne, NX - off - N))], format='lil')
        R[:, ecol] = -evec.reshape(-1, 1); blocks.append(R.tocsr())
    blocks.append(sp.csr_matrix((np.array([1.0, 1.0]), (np.array([0, 0]), np.array([E11, E22]))), shape=(1, NX)))
    b += [0.0] * (3 * ne) + [1.0]; cones.append(clarabel.ZeroConeT(3 * ne + 1))
    blocks.append(sp.csr_matrix((np.array([-1.0, -np.sqrt(2), -1.0]), (np.arange(3), np.array([E11, E12, E22]))), shape=(3, NX)))
    b += [0.0] * 3; cones.append(clarabel.PSDTriangleConeT(2))
    # block maps, one at a time (W blocks, then constraint blocks), as in fr2_herm4.kernels
    specs = [('X', o, m, None) for (o, m) in rel.wcones] + [(name, None, m, Md) for (name, mu, m, Md) in rel.blocks]
    sizes = []; info = []; tot_before = 0; srng = np.random.default_rng(7)
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
        Bt = np.hstack([wp[0] * Np, wm[0] * Nm]); Bb = np.hstack([wp[1] * Np, wm[1] * Nm]); d = dp + dm
        T11 = svec_pair_map(Bt, Bt, svec, svec_pairs) / 2.0; T22 = svec_pair_map(Bb, Bb, svec, svec_pairs) / 2.0
        T12 = svec_pair_map(Bt, Bb, svec, svec_pairs)
        R = np.hstack([T11 @ Mf, T12 @ Mf, T22 @ Mf])            # k_d x 3N, svec of the compressed lifted block
        V = data_span(R, d, unsvec_index(d, svec_pairs), rng=srng); r = V.shape[1]
        tot_before += d
        if r == 0: continue
        Cv = svec_pair_map(V, V, svec, svec_pairs) / 2.0          # svec(V^T D V) from svec(D)
        Rr = Cv @ R; Rr[np.abs(Rr) < 1e-13 * max(1.0, np.abs(Rr).max())] = 0.0
        blocks.append(sp.hstack([sp.csr_matrix(-Rr), sp.csr_matrix((Rr.shape[0], 3))], format='csr'))
        b += [0.0] * Rr.shape[0]; cones.append(clarabel.PSDTriangleConeT(r)); sizes.append(r); info.append((name, dp, dm, r))
        del Mf, R
    A = sp.vstack(blocks, format='csc'); bvec = np.array(b); del blocks
    print(f"   lifted reduced: total face dim {tot_before} -> {sum(sizes)} after data spans; {A.shape[0]} rows, nnz "
          f"{A.nnz / 1e6:.1f}M, {len(sizes)} PSD blocks (largest {max(sizes)}), sum k^2 "
          f"{sum((s * (s + 1) // 2) ** 2 for s in sizes) / 1e6:.1f}M; {time.time() - t0:.0f}s", flush=True)
    c = np.zeros(NX); c[0:N] = rel.m_uu; c[N:2 * N] = 2 * rel.m_uv; c[2 * N:3 * N] = rel.m_vv
    G = (A.T @ A).toarray(); scale = np.sqrt(np.maximum(np.diag(G), 1e-300)); Gs = G / np.outer(scale, scale)
    ev, V = np.linalg.eigh(Gs); null = ev < 1e-12 * ev[-1]; rank = int((~null).sum())
    cs = c / scale; cnull = np.linalg.norm(V[:, null].T @ cs) / np.linalg.norm(cs) if np.any(null) else 0.0
    Lr = V[:, ~null] * np.sqrt(ev[~null])
    _, _, perm = sla.qr(Lr.T, pivoting=True, mode='economic'); J = np.sort(perm[:rank])
    Am = A[:, J]; cobj = c[J]
    print(f"   rank {rank} of {NX} columns; relative |c on null| {cnull:.1e}", flush=True)
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 500
    st.tol_gap_abs = 1e-9; st.tol_gap_rel = 1e-9; st.tol_feas = 1e-9
    try: st.chordal_decomposition_enable = False
    except Exception: pass
    t1 = time.time()
    sol = clarabel.DefaultSolver(sp.csc_matrix((len(J), len(J))), cobj, Am, bvec, cones, st).solve()
    z = np.array(sol.z); r = Am.T @ z + cobj
    row = 3 * ne + 1; zE = z[row:row + 3]; row += 3; mz = np.inf; worst = None; small = 0; eigs = []
    for s, nm in zip(sizes, info):
        kk = s * (s + 1) // 2; e = np.linalg.eigvalsh(unsvec(z[row:row + kk], s)).min(); row += kk; eigs.append(e)
        if e < 1e-6: small += 1
        if e < mz: mz, worst = e, nm
    print(f"   lifted FR ({filt}): status {sol.status}, iters {sol.iterations}, primal {sol.obj_val:+.3e}, "
          f"dual {-(bvec @ z):+.3e}, ||A^T z + c|| {np.abs(r).max():.1e}; min eig reduced lifted dual blocks {mz:.2e} "
          f"(at {worst}); blocks with eig < 1e-6: {small} of {len(sizes)}; dual E-block {np.round(zE, 5).tolist()}; "
          f"solve {time.time() - t1:.0f}s", flush=True)
    np.save('lift_fr2_eigs.npy', np.array(eigs))
