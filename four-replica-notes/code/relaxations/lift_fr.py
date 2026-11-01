"""Facially reduced lifted (all-ratio) certificate at n = 4, in floating point.

Lifted problem (see lift2.py): X = (X11, X12, X22), E = [[e11, e12], [e12, e22]] >= 0, Tr E = 1, equalities lifted,
every PSD block lifted to [[W(X11), W(X12)], [W(X12), W(X22)]], objective m_uu.X11 + 2 m_uv.X12 + m_vv.X22.
A lifted zero is (x (x) w w^T, w w^T) with x a zero at direction w; at t = +1 (projectors) w+ = (1, 1)/sqrt2, at
t = -1 w- = (1, -1)/sqrt2.  Its lifted block is (w w^T) (x) W(x), with range w (x) range W(x).  So the lifted face of
block k is (w+ (x) K+_k) + (w- (x) K-_k), with K+-_k the scalar faces at t = +-1, and the reduced lifted block is
B_k^T [[W11, W12], [W12, W22]] B_k with B_k = [w+ (x) N+_k, w- (x) N-_k].
Zero points are built in parallel worker processes.  Columns outside a maximal independent set of the constraint
matrix (pivoted Cholesky of A^T A) are dropped; the objective is checked to vanish on the null directions.
usage: python lift_fr.py [filter] [n_spinflip] [workers]     filter as in fr2_herm4 (all, COV, ...)"""
import sys, os, time, itertools, numpy as np, scipy.sparse as sp, scipy.linalg as sla, clarabel
os.environ['FR_REAL'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool

_ALG = None
def _zf_worker(args):
    """moments -> zfull tensors for one (u, v) pair (runs in a worker process)"""
    global _ALG
    from alg4 import build_algebras
    from herm4 import FAMS, moments_multi
    from validate import zfull_from_moments
    if _ALG is None: _ALG = build_algebras()
    u, v = args
    return {c: zfull_from_moments(moments_multi([u.reshape(-1, 1) if ch == 'u' else v.reshape(-1, 1) for ch in c], 4), 4, _ALG[0])
            for c in FAMS}

def svec_pair_map(A, B, svec, svec_pairs):
    """matrix of W -> svec(A^T W B + B^T W A) on svec(W) (W symmetric m x m; A, B m x d)"""
    m, d = A.shape; pm = svec_pairs(m); cols = []
    for (i, j) in pm:
        Eij = np.zeros((m, m))
        if i == j: Eij[i, i] = 1.0
        else: Eij[i, j] = Eij[j, i] = 1.0 / np.sqrt(2.0)
        cols.append(svec(A.T @ Eij @ B + B.T @ Eij @ A))
    return np.array(cols).T

def pairs_tp1(rng):
    out = []
    def onb():
        Q, _ = np.linalg.qr(rng.standard_normal((3, 2))); return Q
    for j in range(4):                                   # rotated k1 / k3
        k = 1 if j % 2 == 0 else 3
        p = np.ones(1); q = np.ones(1)
        for s in range(k):
            Q = onb(); p = np.kron(p, Q[:, 0]); q = np.kron(q, Q[:, 1])
        psi = rng.standard_normal(3 ** (4 - k)); psi /= np.linalg.norm(psi)
        u = np.kron(p, psi); v = np.kron(q, psi); th = rng.uniform(0, 2 * np.pi)
        out.append((np.cos(th) * u + np.sin(th) * v, -np.sin(th) * u + np.cos(th) * v))
    for j in range(3):                                   # rotated interpolating
        Qa, Qb, Qc = onb(), onb(), onb(); psi = rng.standard_normal(3); psi /= np.linalg.norm(psi)
        wp = np.kron(np.kron(Qb[:, 0], Qc[:, 0]), psi); wm = np.kron(np.kron(Qb[:, 1], Qc[:, 1]), psi)
        th = rng.uniform(0, np.pi / 2); s = rng.choice([-1, 1])
        phi = np.cos(th) * wp + s * np.sin(th) * wm; chi = np.sin(th) * wp + s * np.cos(th) * wm
        u = np.kron(Qa[:, 0], phi); v = np.kron(Qa[:, 1], chi); th2 = rng.uniform(0, 2 * np.pi)
        out.append((np.cos(th2) * u + np.sin(th2) * v, -np.sin(th2) * u + np.cos(th2) * v))
    return out

def pairs_tm1(rng, nsf):
    from fr_tm1_qubit import spinflip_pair, inner_pair
    out = [inner_pair(rng, True) for j in range(4)]
    out += [spinflip_pair(rng, True) for j in range(nsf)]
    return out

if __name__ == "__main__":
    filt = sys.argv[1] if len(sys.argv) > 1 else 'COV'
    nsf = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    t0 = time.time()
    rng = np.random.default_rng(4242)
    PP = pairs_tp1(rng); PM = pairs_tm1(rng, nsf)
    pool = Pool(nw); job = pool.map_async(_zf_worker, PP + PM)
    from alg4 import build_algebras
    from herm4 import Herm4, ALL_CONS
    from validate import F1_matrix
    from rep4 import unsvec, svec, svec_pairs
    import fr2_herm4 as FR
    keep = FR.FILTERS[filt]
    algs = build_algebras(); F1 = F1_matrix(algs[0])
    rel = Herm4(4, algs, ALL_CONS, verbose=False); N = rel.N
    print(f"n=4 lifted FR ({filt}): build {time.time() - t0:.0f}s ({N} parameters)", flush=True)
    xp = []; xm = []
    for k in (1, 3):
        for j in range(2): xp.append(rel.x_from_zfulls(FR.zero_point_k(4, k, rng, F1))[0])
    for k in (2, 4):
        for j in range(2): xm.append(rel.x_from_zfulls(FR.zero_point_k(4, k, rng, F1))[0])
    zfs = job.get(); pool.close()
    for i, zf in enumerate(zfs):
        x, res = rel.x_from_zfulls(zf)
        (xp if i < len(PP) else xm).append(x)
    print(f"   zeros: {len(xp)} at t=+1, {len(xm)} at t=-1; {time.time() - t0:.0f}s", flush=True)
    for nm, xs, t in (('t=+1', xp, 1.0), ('t=-1', xm, -1.0)):
        worst = max(abs(rel.evaluate(x, 1.0, t)[3]) for x in xs)
        print(f"   {nm}: max |q| {worst:.1e}", flush=True)
    kp = FR.kernels(rel, sum(xp)); km = FR.kernels(rel, sum(xm))
    wp = np.array([1.0, 1.0]) / np.sqrt(2); wm = np.array([1.0, -1.0]) / np.sqrt(2)
    NX = 3 * N + 3; E11, E12, E22 = 3 * N, 3 * N + 1, 3 * N + 2
    Emat = sp.csr_matrix(np.vstack([rel.mom_rows, rel.Q])); evec = np.concatenate([rel.mom_rhs, np.zeros(rel.Q.shape[0])])
    ne = Emat.shape[0]
    blocks = []; b = []; cones = []
    for off, ecol in ((0, E11), (N, E12), (2 * N, E22)):
        R = sp.hstack([sp.csr_matrix((ne, off)), Emat, sp.csr_matrix((ne, NX - off - N))], format='lil')
        R[:, ecol] = -evec.reshape(-1, 1); blocks.append(R.tocsr())
    tr = sp.csr_matrix((np.array([1.0, 1.0]), (np.array([0, 0]), np.array([E11, E22]))), shape=(1, NX))
    blocks.append(tr); b += [0.0] * (3 * ne) + [1.0]; cones.append(clarabel.ZeroConeT(3 * ne + 1))
    blocks.append(sp.csr_matrix((np.array([-1.0, -np.sqrt(2), -1.0]), (np.arange(3), np.array([E11, E12, E22]))), shape=(3, NX)))
    b += [0.0] * 3; cones.append(clarabel.PSDTriangleConeT(2))
    sizes = []; dims = []
    for (name, Mf, Np), (_, _, Nm) in zip(kp, km):
        if name != 'X' and not keep(name): continue
        dp, dm = Np.shape[1], Nm.shape[1]
        if dp + dm == 0: continue
        Bt = np.hstack([wp[0] * Np, wm[0] * Nm]); Bb = np.hstack([wp[1] * Np, wm[1] * Nm])
        T11 = svec_pair_map(Bt, Bt, svec, svec_pairs) / 2.0; T22 = svec_pair_map(Bb, Bb, svec, svec_pairs) / 2.0
        T12 = svec_pair_map(Bt, Bb, svec, svec_pairs)
        d = dp + dm
        R = np.hstack([-(T11 @ Mf), -(T12 @ Mf), -(T22 @ Mf), np.zeros((d * (d + 1) // 2, 3))])
        R[np.abs(R) < 1e-14] = 0.0
        blocks.append(sp.csr_matrix(R)); b += [0.0] * R.shape[0]; cones.append(clarabel.PSDTriangleConeT(d))
        sizes.append(d); dims.append((name, dp, dm))
    del kp, km
    A = sp.vstack(blocks, format='csc'); bvec = np.array(b)
    print(f"   lifted reduced: {A.shape[0]} rows, nnz {A.nnz / 1e6:.1f}M, {len(sizes)} PSD blocks (largest {max(sizes)}), "
          f"sum k^2 {sum((s * (s + 1) // 2) ** 2 for s in sizes) / 1e6:.1f}M; {time.time() - t0:.0f}s", flush=True)
    c = np.zeros(NX); c[0:N] = rel.m_uu; c[N:2 * N] = 2 * rel.m_uv; c[2 * N:3 * N] = rel.m_vv
    G = (A.T @ A).toarray(); scale = np.sqrt(np.maximum(np.diag(G), 1e-300))
    Gs = G / np.outer(scale, scale)
    ev, V = np.linalg.eigh(Gs); null = ev < 1e-12 * ev[-1]; rank = int((~null).sum())
    cnull = np.linalg.norm(V[:, null].T @ (c / scale)) / max(1e-300, np.linalg.norm(c / scale)) if np.any(null) else 0.0
    # maximal independent column set: QR with pivoting on the Cholesky-like factor of Gs
    Lr = V[:, ~null] * np.sqrt(ev[~null])                     # Gs = Lr Lr^T
    _, _, perm = sla.qr(Lr.T, pivoting=True, mode='economic')
    J = np.sort(perm[:rank])
    Am = A[:, J]; cobj = c[J]
    print(f"   rank {rank} of {NX} columns (dropped {NX - rank}); relative |c on null| {cnull:.1e}", flush=True)
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 500
    st.tol_gap_abs = 1e-9; st.tol_gap_rel = 1e-9; st.tol_feas = 1e-9
    try: st.chordal_decomposition_enable = False
    except Exception: pass
    t1 = time.time()
    sol = clarabel.DefaultSolver(sp.csc_matrix((len(J), len(J))), cobj, Am, bvec, cones, st).solve()
    z = np.array(sol.z); r = Am.T @ z + cobj
    row = 3 * ne + 1; zE = z[row:row + 3]; row += 3; mz = np.inf; worst = None; small = 0
    for s, (name, dp, dm) in zip(sizes, dims):
        kk = s * (s + 1) // 2; e = np.linalg.eigvalsh(unsvec(z[row:row + kk], s)).min(); row += kk
        if e < 1e-6: small += 1
        if e < mz: mz, worst = e, (name, dp, dm)
    print(f"   lifted FR ({filt}): status {sol.status}, iters {sol.iterations}, primal {sol.obj_val:+.3e}, "
          f"dual {-(bvec @ z):+.3e}, ||A^T z + c|| {np.abs(r).max():.1e}; min eig reduced lifted dual blocks {mz:.2e} "
          f"(at {worst}); blocks with eig < 1e-6: {small}; dual E-block {np.round(zE, 5).tolist()}; solve {time.time() - t1:.0f}s",
          flush=True)
