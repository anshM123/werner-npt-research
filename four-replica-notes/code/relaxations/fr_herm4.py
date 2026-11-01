"""Dual facial reduction for the two-letter relaxation (herm4) at the endpoint, t = -1 and t = +1.
Known zeros of q_n for C = |u><u| + t |v><v|:
   t = +1:  u = alpha (x) psi, v = alpha' (x) psi   (alpha, alpha' orthonormal on one site, psi pure on the others)
   t = -1:  u = a (x) b (x) psi, v = a' (x) b' (x) psi   (a ⊥ a', b ⊥ b' on two sites, psi pure on the others)
(the local search at n = 2, 3 lands on the second family).  Their two-letter moments X_c factorise over the site groups,
so they are computed from small moment tensors.  Every certificate of value 0 must vanish on the ranges of all PSD blocks
at these points; we compress every block to the kernel at the sum of K sampled zeros and re-solve.
usage: python fr_herm4.py n t K [L]   (L = number of extra zeros taken from local-search minimisers)"""
import sys, time, itertools, numpy as np, scipy.sparse as sp, clarabel
from alg4 import build_algebras
from herm4 import Herm4, FAMS, moments_multi, ALL_CONS
from validate import F1_matrix
from rep4 import unsvec, svec, svec_pairs

def hat_unsym(Us, n, F1):
    """unsymmetrised hat-coefficient tensor (23,)*n of twirl((x)_r U_r U_r^dagger)"""
    Z = moments_multi(Us, n)
    for s in range(n): Z = np.moveaxis(np.tensordot(F1, Z, axes=([1], [s])), 0, s)
    return Z.real

def rand_unit(rng, D):
    x = rng.standard_normal(D) + 1j * rng.standard_normal(D); return (x / np.linalg.norm(x)).reshape(D, 1)

def zero_point(n, t, rng, F1):
    """dict family -> site-symmetrised hat tensor for one random zero"""
    zf = {}
    if t > 0:
        Q, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
        al, al2 = Q[:, :1], Q[:, 1:2]; psi = rand_unit(rng, 3 ** (n - 1))
        for c in FAMS:
            z1 = hat_unsym([al if ch == 'u' else al2 for ch in c], 1, F1)
            zr = hat_unsym([psi] * 4, n - 1, F1)
            zf[c] = np.multiply.outer(z1, zr)
    else:
        Qa, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
        Qb, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
        psi = rand_unit(rng, 3 ** (n - 2)) if n > 2 else None
        for c in FAMS:
            za = hat_unsym([Qa[:, :1] if ch == 'u' else Qa[:, 1:2] for ch in c], 1, F1)
            zb = hat_unsym([Qb[:, :1] if ch == 'u' else Qb[:, 1:2] for ch in c], 1, F1)
            T = np.multiply.outer(za, zb)
            if psi is not None: T = np.multiply.outer(T, hat_unsym([psi] * 4, n - 2, F1))
            zf[c] = T
    perms = list(itertools.permutations(range(n)))
    return {c: sum(T.transpose(p) for p in perms) / len(perms) for c, T in zf.items()}

def svec_compress(Nk):
    m, k = Nk.shape; pm = svec_pairs(m)
    Cm = np.zeros((k * (k + 1) // 2, len(pm)))
    for b_, (i, j) in enumerate(pm):
        E = np.zeros((m, m))
        if i == j: E[i, i] = 1.0
        else: E[i, j] = E[j, i] = 1.0 / np.sqrt(2.0)
        Cm[:, b_] = svec(Nk.T @ E @ Nk)
    return Cm

if __name__ == "__main__":
    n = int(sys.argv[1]); t = float(sys.argv[2]); K = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    algs = build_algebras(); F1 = F1_matrix(algs[0])
    t0 = time.time()
    rel = Herm4(n, algs, ALL_CONS, verbose=False)
    print(f"build {time.time() - t0:.0f}s ({rel.N} parameters)", flush=True)
    rng = np.random.default_rng(31)
    xs = []
    for k in range(K):
        zf = zero_point(n, t, rng, F1)
        x, res = rel.x_from_zfulls(zf)
        mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
        print(f"   zero {k}: reconstruction {res:.1e}, q = {obj:.1e}, moment res {mom:.1e}, identity res {eqr:.1e}, "
              f"min eig {min(mins.values()):.1e}", flush=True)
        xs.append(x)
    L = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    if L:
        from local_herm import run
        from validate import zfull_from_moments
        seed = 100; got = 0
        while got < L and seed < 100 + 6 * L:
            best, U = run(n, t, restarts=1, iters=4000, seed=seed); seed += 1
            if best > 1e-11: continue
            u, v = U[:, :1], U[:, 1:2]
            zf = {c: zfull_from_moments(moments_multi([u if ch == 'u' else v for ch in c], n), n, algs[0]) for c in FAMS}
            x, res = rel.x_from_zfulls(zf)
            mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
            print(f"   local-search zero (seed {seed - 1}, q {best:.1e}): reconstruction {res:.1e}, q = {obj:.1e}, "
                  f"min eig {min(mins.values()):.1e}; {time.time() - t0:.0f}s", flush=True)
            xs.append(x); got += 1
    xstar = sum(xs)
    # full row matrices of all PSD constraints (W blocks first), compressed to the kernels at xstar
    N = rel.N; E = np.vstack([rel.mom_rows, rel.Q]); b = list(rel.mom_rhs) + [0.0] * rel.Q.shape[0]
    rows = [E]; cones = [clarabel.ZeroConeT(E.shape[0])]; psd = []; nn = []; kept = 0; dropped = 0; dims = [0, 0]
    def handle(Mfull, m):
        global kept, dropped
        X = unsvec(Mfull @ xstar, m); ev, V = np.linalg.eigh(X)
        Nk = V[:, ev < 1e-9 * max(1.0, ev.max())]; dims[0] += m; dims[1] += Nk.shape[1]
        if Nk.shape[1] == 0: dropped += 1; return
        kept += 1; Cm = svec_compress(Nk) @ Mfull
        if Nk.shape[1] == 1: nn.append(-Cm[0])
        else: psd.append((Nk.shape[1], -Cm))
    for (o, m) in rel.wcones:
        k = m * (m + 1) // 2; Mf = np.zeros((k, N)); Mf[np.arange(k), o + np.arange(k)] = 1.0; handle(Mf, m)
    for (name, mu, m, Md) in rel.blocks:
        k = m * (m + 1) // 2; Mf = np.zeros((k, N))
        for c, Mc in Md.items(): Mf[:, rel.fam[c]['start']:rel.fam[c]['stop']] = Mc
        handle(Mf, m)
    if nn: rows.append(np.array(nn)); b += [0.0] * len(nn); cones.append(clarabel.NonnegativeConeT(len(nn)))
    for k, Cm in psd: rows.append(Cm); b += [0.0] * Cm.shape[0]; cones.append(clarabel.PSDTriangleConeT(k))
    Ad = np.vstack(rows); bvec = np.array(b)
    print(f"   blocks kept {kept}, dropped {dropped}; block dimension {dims[0]} -> {dims[1]}; reduced A {Ad.shape}; "
          f"{time.time() - t0:.0f}s", flush=True)
    # row space of Ad from the Gram matrix Ad^T Ad (N x N); null directions are checked directly on Ad
    Gm = np.zeros((N, N))
    for i0 in range(0, Ad.shape[0], 20000): Gm += Ad[i0:i0 + 20000].T @ Ad[i0:i0 + 20000]
    ev_, V_ = np.linalg.eigh(Gm)
    null = ev_ < 1e-13 * ev_[-1]
    if np.any(null):
        chk = np.linalg.norm(Ad @ V_[:, null], axis=0).max() / np.sqrt(ev_[-1])
        print(f"   null directions {int(null.sum())}, max |Ad v| / |Ad| = {chk:.1e}, smallest kept eigenvalue ratio "
              f"{ev_[~null].min() / ev_[-1]:.1e}", flush=True)
    rk = int((~null).sum()); Rb = V_[:, ~null].T
    c = 1.0 * rel.m_uu + t * t * rel.m_vv + 2 * t * rel.m_uv
    Nul = np.eye(N) - Rb.T @ Rb
    print(f"   rank of reduced A {rk} / {N}; |objective on null(A)| = {np.linalg.norm(Nul @ c):.1e}", flush=True)
    Amat = sp.csc_matrix(Ad @ Rb.T); cobj = Rb @ c
    for tol in (1e-9, 1e-11):
        st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 500
        st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
        t1 = time.time()
        sol = clarabel.DefaultSolver(sp.csc_matrix((rk, rk)), cobj, Amat, bvec, cones, st).solve()
        z = np.array(sol.z); r = Amat.T @ z + cobj
        row = E.shape[0]; mz = np.inf
        if nn: mz = min(mz, z[row:row + len(nn)].min()); row += len(nn)
        for k, Cm in psd:
            kk = Cm.shape[0]; mz = min(mz, np.linalg.eigvalsh(unsvec(z[row:row + kk], k)).min()); row += kk
        print(f"   reduced problem (tol {tol:g}): status {sol.status}, iters {sol.iterations}, primal q {sol.obj_val:.4e}, "
              f"dual q {-(bvec @ z):.4e}, ||A^T z + c|| {np.abs(r).max():.1e}, min eig reduced dual blocks {mz:.1e}; "
              f"{time.time() - t1:.0f}s", flush=True)
