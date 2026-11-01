"""n = 4 checks on the swap-reduced two-letter relaxation (herm4):
 (b) true-point validation on the n = 4 build: random orthonormal pairs (complex and real), every block PSD, moments and
     identities exact, objective = q_4(C) computed directly;
 (a) dual lower bounds of the grid solutions saved by grid_herm4.py: for each t, with the saved (x, z, s):
     consistency of the rebuilt data (|b - A x - s|), dual objective -b.z, residual r = A^T z + c, ||r||_inf,
     ||x||_1, -b.z - ||r||_inf ||x||_1, and the smallest eigenvalue of the dual blocks.
usage: python check_n4.py [ntrials]"""
import sys, time, glob, numpy as np
from alg4 import build_algebras
from herm4 import Herm4, FAMS, moments_multi, ALL_CONS
from validate import zfull_from_moments, qn_direct
from rep4 import unsvec

n = 4; ntrials = int(sys.argv[1]) if len(sys.argv) > 1 else 2
algs = build_algebras()
t0 = time.time()
rel = Herm4(n, algs, ALL_CONS, verbose=False)
print(f"build {time.time() - t0:.0f}s: {rel.N} parameters, identity rank {rel.Q.shape[0]}, {len(rel.blocks)} PSD blocks", flush=True)
# ---------------- (b) true points ----------------
rng = np.random.default_rng(404)
for trial in range(ntrials):
    D = 3 ** n
    X = rng.standard_normal((D, 2)) + (1j * rng.standard_normal((D, 2)) if trial % 2 == 0 else 0)
    Q, _ = np.linalg.qr(X); u = Q[:, :1]; v = Q[:, 1:2]
    t1 = time.time()
    zf = {c: zfull_from_moments(moments_multi([u if ch == 'u' else v for ch in c], n), n, algs[0]) for c in FAMS}
    x, res = rel.x_from_zfulls(zf)
    for (cu, cv) in [(1.0, -1.0), (1.0, 1.0), (1.0, -0.25)]:
        mins, mom, eqr, obj = rel.evaluate(x, cu, cv)
        C = cu * (u @ u.conj().T) + cv * (v @ v.conj().T)
        print(f"(b) n=4 trial {trial} ({'complex' if trial % 2 == 0 else 'real'}) C=({cu},{cv}): reconstruction {res:.1e}, "
              f"moment res {mom:.1e}, identity res {eqr:.1e}, objective {obj:.12f} vs q_4(C) {qn_direct(C, n):.12f}; "
              f"min eig {', '.join(f'{k}:{v_:.1e}' for k, v_ in mins.items())}; {time.time() - t1:.0f}s", flush=True)
# ---------------- (a) dual bounds of the saved grid solutions ----------------
A, b, cones = rel.assemble(); rel.Aasm = (A, b, cones)
def walk_min(z):
    row = rel.mom_rows.shape[0] + rel.Q.shape[0]; mz = np.inf
    nn = sum(1 for (o, m) in rel.wcones if m == 1) + sum(1 for blk in rel.blocks if blk[2] == 1)
    if nn: mz = min(mz, z[row:row + nn].min()); row += nn
    for m in [m for (o, m) in rel.wcones if m > 1] + [blk[2] for blk in rel.blocks if blk[2] > 1]:
        k = m * (m + 1) // 2; mz = min(mz, np.linalg.eigvalsh(unsvec(z[row:row + k], m)).min()); row += k
    assert row == A.shape[0]
    return mz
for fn in sorted(glob.glob('herm4_n4_all_t*.npz')):
    dat = np.load(fn); x = dat['x']; z = dat['z']; s = dat['s']
    t = float(fn.split('_t')[1][:-4])
    c = 1.0 * rel.m_uu + t * t * rel.m_vv + 2 * t * rel.m_uv
    r = A.T @ z + c
    cons_err = np.abs(b - A @ x - s).max()
    dual = -(b @ z); prim = c @ x
    true = (1 - abs(t)) ** 2 / 16
    print(f"(a) t={t:+.2f}: primal c.x {prim:.8e}, dual -b.z {dual:.8e} (true min {true:.8e}; dual - true {dual - true:+.1e}), "
          f"||A^T z + c|| {np.abs(r).max():.1e}, ||x||_1 {np.abs(x).sum():.2f}, corrected {dual - np.abs(r).max() * np.abs(x).sum():.2e}, "
          f"min eig dual blocks {walk_min(z):.1e}, rebuild consistency |b - A x - s| {cons_err:.1e}", flush=True)
