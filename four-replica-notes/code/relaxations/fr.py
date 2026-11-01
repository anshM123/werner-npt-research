"""Dual facial reduction test at the endpoint.
Known optimal points (q_n = 0): Z* = site-symmetrised twirl of (P_1 (x) |0><0|^(x)(n-1))^(x)4, P_1 a rank-2 projector on
one site.  If the relaxation is exact, some optimal dual certificate vanishes on the range of every PSD block at Z*
(complementary slackness).  We therefore replace each PSD constraint X(w) >= 0 by N^T X(w) N >= 0, N = orthonormal
basis of ker X(w*), and drop blocks with trivial kernel.  The reduced primal has a larger feasible set; its optimum is 0
iff a certificate of value 0 exists that is supported on these kernels (and it is unbounded/negative otherwise).
usage: python fr.py n cons [w]"""
import sys, time, itertools, numpy as np, scipy.sparse as sp, clarabel
from alg4 import *
from rep4 import Relax, unsvec, svec, svec_pairs

def site_hat(Ps, A):
    P4 = Ps
    for _ in range(3): P4 = np.kron(P4, Ps)
    out = []
    for lam in range(A.nb):
        f = A.f[lam]
        for i in range(f):
            for j in range(f):
                out.append(np.trace(P4 @ A.unit(lam, j, i)).real)
    return np.array(out)

def zstar_full(n, A, K=0, seed=5):
    """site-symmetrised sum of twirls of P_1 (x) |psi><psi| (psi pure on n-1 sites): psi = |0..0> plus K random psi"""
    from validate import moments, F1_matrix
    P1 = np.diag([1.0, 1.0, 0.0]); P0 = np.diag([1.0, 0.0, 0.0])
    z1, z0 = site_hat(P1, A), site_hat(P0, A)
    T0 = z0
    for _ in range(n - 2): T0 = np.multiply.outer(T0, z0)
    tails = [T0]
    rng = np.random.default_rng(seed); F1 = F1_matrix(A)
    for k in range(K):
        m = n - 1; v = rng.standard_normal(3 ** m) + 1j * rng.standard_normal(3 ** m); v /= np.linalg.norm(v)
        mom = moments(np.outer(v, v.conj()), m)
        Z = mom
        for s in range(m):
            Z = np.moveaxis(np.tensordot(F1, Z, axes=([1], [s])), 0, s)
        tails.append(Z.real)
    T = sum(np.multiply.outer(z1, t) for t in tails)
    perms = list(itertools.permutations(range(n)))
    return sum(T.transpose(p) for p in perms) / len(perms)

def svec_compress(N):
    """matrix C with svec(N^T X N) = C svec(X) for symmetric X (m x m), N (m x k)"""
    m, k = N.shape
    pm = svec_pairs(m); pk = svec_pairs(k)
    C = np.zeros((len(pk), len(pm)))
    for b, (i, j) in enumerate(pm):
        E = np.zeros((m, m))
        if i == j: E[i, i] = 1.0
        else: E[i, j] = E[j, i] = 1.0 / np.sqrt(2.0)
        C[:, b] = svec(N.T @ E @ N)
    return C

if __name__ == "__main__":
    n = int(sys.argv[1]); cons = sys.argv[2].split(','); w = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    RANK = 2.0
    algs = build_algebras(); A = algs[0]
    rel = Relax(n, RANK, algs, cons, verbose=False, w=w)
    K = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    Zf = zstar_full(n, A, K)
    ws, res = rel.w_from_zfull(Zf)
    mins, mres, eres, objv = rel.evaluate(ws)
    print(f"n={n} cons={cons} w={w}: Z* reconstruction residual {res:.1e}, moment res {mres:.1e}, objective q(Z*) = {objv / 4:.3e}, "
          f"min eigs {', '.join(f'{k}:{v:.1e}' for k, v in mins.items())}", flush=True)
    # collect blocks: (row-matrix M mapping w -> svec(X), size m)
    blocks = []
    for (o, m) in rel.wcones:
        k = m * (m + 1) // 2; M = np.zeros((k, rel.N)); M[np.arange(k), o + np.arange(k)] = 1.0
        blocks.append(('W', m, M))
    for (cn, mu, c, M) in rel.blocks: blocks.append((cn, c['m'], M))
    Arows = [np.array(rel.mom_rows)]; b = list(rel.mom_rhs); cones = [clarabel.ZeroConeT(len(rel.mom_rhs))]
    if rel.eq_rows.shape[0]:
        Arows.append(rel.eq_rows); b += [0.0] * rel.eq_rows.shape[0]; cones[0] = clarabel.ZeroConeT(len(b))
    stats = {'kept': 0, 'dropped': 0, 'dim_before': 0, 'dim_after': 0}
    nn_rows = []; psd = []
    for (name, m, M) in blocks:
        X = unsvec(M @ ws, m)
        ev, V = np.linalg.eigh(X)
        Nk = V[:, ev < 1e-9 * max(1.0, ev.max())]
        stats['dim_before'] += m
        if Nk.shape[1] == 0:
            stats['dropped'] += 1; continue
        stats['kept'] += 1; stats['dim_after'] += Nk.shape[1]
        C = svec_compress(Nk) @ M
        if Nk.shape[1] == 1: nn_rows.append(-C[0])
        else: psd.append((Nk.shape[1], -C))
    if nn_rows:
        Arows.append(np.array(nn_rows)); b += [0.0] * len(nn_rows); cones.append(clarabel.NonnegativeConeT(len(nn_rows)))
    for k, C in psd:
        Arows.append(C); b += [0.0] * C.shape[0]; cones.append(clarabel.PSDTriangleConeT(k))
    print(f"   blocks kept {stats['kept']} dropped {stats['dropped']}; total block dim {stats['dim_before']} -> {stats['dim_after']}", flush=True)
    Amat = sp.csc_matrix(np.vstack(Arows)); bvec = np.array(b)
    Ad = np.vstack(Arows)
    u_, s_, vt_ = np.linalg.svd(Ad, full_matrices=True)
    rk = int(np.sum(s_ > 1e-9 * s_[0])); Nul = vt_[rk:]
    print(f"   reduced A: {Ad.shape}, rank {rk}; |projection of objective on null(A)| = "
          f"{np.linalg.norm(Nul @ rel.obj) if Nul.shape[0] else 0.0:.2e}", flush=True)
    if Nul.shape[0]:
        # restrict to the row space of A (directions in null(A) are unconstrained)
        Rb = vt_[:rk]                      # w = Rb^T y
        Amat = sp.csc_matrix(Ad @ Rb.T); cobj = Rb @ rel.obj
    else:
        Rb = np.eye(rel.N); cobj = rel.obj.copy()
    for tol in (1e-9, 1e-11):
        st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 1000
        st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
        t0 = time.time()
        sol = clarabel.DefaultSolver(sp.csc_matrix((Amat.shape[1], Amat.shape[1])), cobj.copy(), Amat, bvec, cones, st).solve()
        x = Rb.T @ np.array(sol.x); z = np.array(sol.z); r = Amat.T @ z + cobj
        print(f"   reduced problem (tol {tol:g}): status {sol.status}, iters {sol.iterations}, primal q = {rel.obj @ x / 4:.4e}, "
              f"dual q = {-(bvec @ z) / 4:.4e}, ||A^T z + c|| = {np.abs(r).max():.1e}; {time.time() - t0:.1f}s", flush=True)
        # min eigenvalue of dual blocks (strict feasibility of the reduced certificate)
        row = len(b) - sum(C.shape[0] for _, C in psd) - len(nn_rows)
        mz = np.inf
        if nn_rows: mz = min(mz, z[row:row + len(nn_rows)].min()); row += len(nn_rows)
        for k, C in psd:
            kk = C.shape[0]; mz = min(mz, np.linalg.eigvalsh(unsvec(z[row:row + kk], k)).min()); row += kk
        print(f"      min eigenvalue over reduced dual blocks: {mz:.2e}", flush=True)
        np.savez(f"fr_n{n}_{'_'.join(cons)}_K{K}_tol{tol:g}.npz", x=x, y=np.array(sol.x), z=z, Rb=Rb, wstar=ws, obj=rel.obj)
