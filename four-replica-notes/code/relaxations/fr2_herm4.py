"""Facial reduction of the two-letter relaxation (herm4) at the endpoint with structured zero families and
local-search zeros, plus a constraint-family ablation on the reduced problem.

Zero families of q_n for C = |u><u| + t|v><v| (checked directly): u, v orthogonal products on a set S of k sites,
sharing a pure state psi (entangled when n - k >= 2) on the other sites:  u = (x)_{s in S} a_s (x) psi,
v = (x)_{s in S} a'_s (x) psi with a_s ⊥ a'_s.  For such C,
   q_n = (1 + t^2)[(1/2)^k - (-1/2)^k] q_{n-k}(psi) + (-1/2)^k (1 + t)^2 q_{n-k}(psi) ... = 0
exactly when k is odd and t = +1, or k is even and t = -1.
usage: python fr2_herm4.py n t ks K L variants
   ks: comma list of k values; K: samples per k; L: local-search zeros; variants: ';'-list of block filters
   (all, noPT2, COV+PT1, COV, PT1, PT, noCOV, noLE2, COV+LE2, PT1+LE2)"""
import sys, time, itertools, numpy as np, scipy.sparse as sp, clarabel
from alg4 import build_algebras
from herm4 import Herm4, FAMS, moments_multi, ALL_CONS
from validate import F1_matrix, zfull_from_moments
from rep4 import unsvec, svec, svec_pairs

def hat_unsym(Us, n, F1):
    Z = moments_multi(Us, n)
    for s in range(n): Z = np.moveaxis(np.tensordot(F1, Z, axes=([1], [s])), 0, s)
    return Z.real

def rand_unit(rng, D):
    x = rng.standard_normal(D) + 1j * rng.standard_normal(D); return (x / np.linalg.norm(x)).reshape(D, 1)

def zero_point_k(n, k, rng, F1):
    """orthogonal products on k sites, shared random psi on the other n - k sites; site-symmetrised hat tensors"""
    pairs = []
    for s in range(k):
        Q, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2))); pairs.append((Q[:, :1], Q[:, 1:2]))
    psi = rand_unit(rng, 3 ** (n - k)) if n - k > 0 else None
    zf = {}
    for c in FAMS:
        T = None
        for (a, a2) in pairs:
            z = hat_unsym([a if ch == 'u' else a2 for ch in c], 1, F1)
            T = z if T is None else np.multiply.outer(T, z)
        if psi is not None:
            z = hat_unsym([psi] * 4, n - k, F1)
            T = z if T is None else np.multiply.outer(T, z)
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

FILTERS = {'all': lambda nm: True, 'noPT2': lambda nm: not nm.startswith('PT2'),
           'COV+PT1': lambda nm: nm.startswith('COV') or nm.startswith('PT1'), 'COV': lambda nm: nm.startswith('COV'),
           'PT1': lambda nm: nm.startswith('PT1'), 'PT': lambda nm: nm.startswith('PT'),
           'noCOV': lambda nm: not nm.startswith('COV'), 'noLE2': lambda nm: not nm.startswith('LE2'),
           'COV+LE2': lambda nm: nm.startswith('COV') or nm.startswith('LE2'),
           'PT1+LE2': lambda nm: nm.startswith('PT1') or nm.startswith('LE2'), 'none': lambda nm: False}

def kernels(rel, xstar):
    """per PSD block (W blocks first, then constraint blocks): (name, full row matrix, kernel basis at xstar)"""
    N = rel.N; out = []
    for (o, m) in rel.wcones:
        k = m * (m + 1) // 2; Mf = np.zeros((k, N)); Mf[np.arange(k), o + np.arange(k)] = 1.0
        out.append(('X', Mf, m))
    for (name, mu, m, Md) in rel.blocks:
        k = m * (m + 1) // 2; Mf = np.zeros((k, N))
        for c, Mc in Md.items(): Mf[:, rel.fam[c]['start']:rel.fam[c]['stop']] = Mc
        out.append((name, Mf, m))
    res = []
    for (name, Mf, m) in out:
        X = unsvec(Mf @ xstar, m); ev, V = np.linalg.eigh(X)
        Nk = V[:, ev < 1e-9 * max(1.0, ev.max())]
        res.append((name, Mf, Nk))
    return res

def solve_reduced(rel, t, kers, keep, tol=1e-9, label=''):
    N = rel.N; E = np.vstack([rel.mom_rows, rel.Q]); b = list(rel.mom_rhs) + [0.0] * rel.Q.shape[0]
    rows = [E]; cones = [clarabel.ZeroConeT(E.shape[0])]; nn = []; psd = []; dims = [0, 0]
    for (name, Mf, Nk) in kers:
        if name != 'X' and not keep(name): continue
        dims[0] += Mf.shape[0]
        if Nk.shape[1] == 0: continue
        Cm = svec_compress(Nk) @ Mf; dims[1] += Cm.shape[0]
        if Nk.shape[1] == 1: nn.append(-Cm[0])
        else: psd.append((Nk.shape[1], -Cm))
    if nn: rows.append(np.array(nn)); b += [0.0] * len(nn); cones.append(clarabel.NonnegativeConeT(len(nn)))
    for k, Cm in psd: rows.append(Cm); b += [0.0] * Cm.shape[0]; cones.append(clarabel.PSDTriangleConeT(k))
    Ad = np.vstack(rows); bvec = np.array(b)
    Gm = np.zeros((N, N))
    for i0 in range(0, Ad.shape[0], 20000): Gm += Ad[i0:i0 + 20000].T @ Ad[i0:i0 + 20000]
    ev_, V_ = np.linalg.eigh(Gm); null = ev_ < 1e-13 * ev_[-1]
    c = 1.0 * rel.m_uu + t * t * rel.m_vv + 2 * t * rel.m_uv
    Rb = V_[:, ~null].T; cnull = np.linalg.norm(V_[:, null].T @ c) if np.any(null) else 0.0
    Amat = sp.csc_matrix(Ad @ Rb.T); cobj = Rb @ c
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 500
    st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
    t1 = time.time()
    sol = clarabel.DefaultSolver(sp.csc_matrix((Rb.shape[0], Rb.shape[0])), cobj, Amat, bvec, cones, st).solve()
    z = np.array(sol.z); r = Amat.T @ z + cobj
    row = E.shape[0]; mz = np.inf
    if nn: mz = min(mz, z[row:row + len(nn)].min()); row += len(nn)
    for k, Cm in psd:
        kk = Cm.shape[0]; mz = min(mz, np.linalg.eigvalsh(unsvec(z[row:row + kk], k)).min()); row += kk
    print(f"   [{label}] reduced rows {dims[0]} -> {dims[1]}, null dirs {int(null.sum())} (|c on null| {cnull:.1e}); "
          f"status {sol.status}, iters {sol.iterations}, primal q {sol.obj_val:.4e}, dual q {-(bvec @ z):.4e}, "
          f"||A^T z + c|| {np.abs(r).max():.1e}, min eig reduced dual blocks {mz:.1e}; {time.time() - t1:.0f}s", flush=True)
    return sol

if __name__ == "__main__":
    n = int(sys.argv[1]); t = float(sys.argv[2]); ks = [int(a) for a in sys.argv[3].split(',')]
    K = int(sys.argv[4]); L = int(sys.argv[5]); variants = sys.argv[6].split(';') if len(sys.argv) > 6 else ['all']
    algs = build_algebras(); F1 = F1_matrix(algs[0])
    t0 = time.time()
    rel = Herm4(n, algs, ALL_CONS, verbose=False)
    print(f"n={n} t={t:+.0f}: build {time.time() - t0:.0f}s ({rel.N} parameters, {len(rel.blocks)} PSD blocks)", flush=True)
    rng = np.random.default_rng(77)
    xs = []
    for k in ks:
        for j in range(K):
            zf = zero_point_k(n, k, rng, F1); x, res = rel.x_from_zfulls(zf)
            mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
            print(f"   zero k={k} #{j}: reconstruction {res:.1e}, q {obj:.1e}, moment res {mom:.1e}, identity res {eqr:.1e}, "
                  f"min eig {min(mins.values()):.1e}", flush=True)
            xs.append(x)
    # optional rotated-basis zeros (t = +1 only): computed in-process (R per k) or loaded from saved x vectors
    import glob, os
    R = int(os.environ.get('FR_ROT', '0')); rot_glob = os.environ.get('FR_ROT_FILES', '')
    if t > 0 and R:
        for k in ks:
            for j in range(R):
                u = np.ones(1, complex); v = np.ones(1, complex)
                for s in range(k):
                    Qp, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
                    u = np.kron(u, Qp[:, 0]); v = np.kron(v, Qp[:, 1])
                if n - k > 0:
                    psi = rand_unit(rng, 3 ** (n - k)).ravel(); u = np.kron(u, psi); v = np.kron(v, psi)
                Rm, _ = np.linalg.qr(rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2)))
                Um = np.stack([u, v], axis=1) @ Rm
                zf = {c: zfull_from_moments(moments_multi([Um[:, :1] if ch == 'u' else Um[:, 1:2] for ch in c], n), n, algs[0]) for c in FAMS}
                x, res = rel.x_from_zfulls(zf); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
                print(f"   rotated zero k={k} #{j}: reconstruction {res:.1e}, q {obj:.1e}, min eig {min(mins.values()):.1e}", flush=True)
                xs.append(x)
    # interpolating family at t = +1: u = a (x) phi, v = a' (x) chi with phi = cos th w+ + e^{i al} sin th w-,
    # chi = sin th w+ + e^{i al} cos th w-, (w+, w-) a t = -1 zero on the other n - 1 sites (orthogonal products on an
    # even number k' of sites sharing psi on the rest); then a random rotation of the (u, v) basis
    Fn = int(os.environ.get('FR_FNEW', '0'))
    if t > 0 and Fn:
        for j in range(Fn):
            kp = 2
            wp = np.ones(1, complex); wm = np.ones(1, complex)
            for s in range(kp):
                Qp, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
                wp = np.kron(wp, Qp[:, 0]); wm = np.kron(wm, Qp[:, 1])
            if n - 1 - kp > 0:
                psi = rand_unit(rng, 3 ** (n - 1 - kp)).ravel(); wp = np.kron(wp, psi); wm = np.kron(wm, psi)
            th = rng.uniform(0, np.pi / 2); al = rng.uniform(0, 2 * np.pi)
            phi = np.cos(th) * wp + np.exp(1j * al) * np.sin(th) * wm; chi = np.sin(th) * wp + np.exp(1j * al) * np.cos(th) * wm
            Qa, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
            u = np.kron(Qa[:, 0], phi); v = np.kron(Qa[:, 1], chi)
            Rm, _ = np.linalg.qr(rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2)))
            Um = np.stack([u, v], axis=1) @ Rm
            zf = {c: zfull_from_moments(moments_multi([Um[:, :1] if ch == 'u' else Um[:, 1:2] for ch in c], n), n, algs[0]) for c in FAMS}
            x, res = rel.x_from_zfulls(zf); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
            print(f"   interpolating zero #{j} (theta {th:.2f}): reconstruction {res:.1e}, q {obj:.1e}, min eig {min(mins.values()):.1e}", flush=True)
            xs.append(x)
    # t = -1 family: u = a (x) phi, v = a' (x) chi, (phi, chi) an arbitrary orthonormal basis (random 2x2 unitary)
    # of the range of a t = +1 zero projector on the other n - 1 sites (orthogonal products on 3 sites, or the
    # interpolating family); then q_n(|u><u| - |v><v|) = q_{n-1}(phi phi^+ + chi chi^+) / 2 = 0
    Tm = int(os.environ.get('FR_TM1NEW', '0'))
    if t < 0 and Tm:
        def onb_pairs(m_sites, kind):
            if kind == 'orth':
                p = np.ones(1, complex); q_ = np.ones(1, complex)
                kk = 3 if m_sites >= 3 else 1
                for s in range(kk):
                    Qp, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
                    p = np.kron(p, Qp[:, 0]); q_ = np.kron(q_, Qp[:, 1])
                if m_sites - kk > 0:
                    psi = rand_unit(rng, 3 ** (m_sites - kk)).ravel(); p = np.kron(p, psi); q_ = np.kron(q_, psi)
                return p, q_
            # interpolating: a (x) (cos th w+ + e^{i al} sin th w-),  a' (x) (sin th w+ + e^{i al} cos th w-)
            wp = np.ones(1, complex); wm = np.ones(1, complex)
            for s in range(2):
                Qp, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
                wp = np.kron(wp, Qp[:, 0]); wm = np.kron(wm, Qp[:, 1])
            if m_sites - 3 > 0:
                psi = rand_unit(rng, 3 ** (m_sites - 3)).ravel(); wp = np.kron(wp, psi); wm = np.kron(wm, psi)
            th = rng.uniform(0, np.pi / 2); al = rng.uniform(0, 2 * np.pi)
            Qa, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
            return (np.kron(Qa[:, 0], np.cos(th) * wp + np.exp(1j * al) * np.sin(th) * wm),
                    np.kron(Qa[:, 1], np.sin(th) * wp + np.exp(1j * al) * np.cos(th) * wm))
        for j in range(Tm):
            kind = 'orth' if (j % 2 == 0 or n - 1 < 3) else 'interp'
            p, q_ = onb_pairs(n - 1, kind)
            Rm, _ = np.linalg.qr(rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2)))
            PQ = np.stack([p, q_], axis=1) @ Rm
            Qa, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
            zf = {}
            for c in FAMS:
                z1 = hat_unsym([Qa[:, :1] if ch == 'u' else Qa[:, 1:2] for ch in c], 1, F1)
                z2 = hat_unsym([PQ[:, :1] if ch == 'u' else PQ[:, 1:2] for ch in c], n - 1, F1)
                zf[c] = np.multiply.outer(z1, z2)
            perms = list(itertools.permutations(range(n)))
            zf = {c: sum(T.transpose(pp) for pp in perms) / len(perms) for c, T in zf.items()}
            x, res = rel.x_from_zfulls(zf); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
            print(f"   t=-1 inner-{kind} zero #{j}: reconstruction {res:.1e}, q {obj:.1e}, min eig {min(mins.values()):.1e}", flush=True)
            xs.append(x)
    if rot_glob:
        for fn in sorted(glob.glob(rot_glob)):
            x = np.load(fn); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
            print(f"   loaded zero {fn}: q {obj:.1e}, moment res {mom:.1e}, identity res {eqr:.1e}, min eig {min(mins.values()):.1e}", flush=True)
            xs.append(x)
    xstar = sum(xs)
    kers = kernels(rel, xstar)
    print(f"   zeros used: {len(xs)}; {time.time() - t0:.0f}s", flush=True)
    if not os.environ.get('FR_SKIP_FIRST'):
        for v in variants: solve_reduced(rel, t, kers, FILTERS[v], label=f'structured, {v}')
    if L:
        from local_herm import run
        seed = 200; got = 0
        while got < L and seed < 200 + 8 * L:
            best, U = run(n, t, restarts=1, iters=4000, seed=seed); seed += 1
            if best > 1e-11: continue
            u, v = U[:, :1], U[:, 1:2]
            zf = {c: zfull_from_moments(moments_multi([u if ch == 'u' else v for ch in c], n), n, algs[0]) for c in FAMS}
            x, res = rel.x_from_zfulls(zf); mins, mom, eqr, obj = rel.evaluate(x, 1.0, t)
            print(f"   local-search zero (seed {seed - 1}, q {best:.1e}): reconstruction {res:.1e}, q {obj:.1e}, "
                  f"min eig {min(mins.values()):.1e}; {time.time() - t0:.0f}s", flush=True)
            xs.append(x); got += 1
        kers = kernels(rel, sum(xs))
        for v in variants: solve_reduced(rel, t, kers, FILTERS[v], label=f'structured + {got} local-search, {v}')
