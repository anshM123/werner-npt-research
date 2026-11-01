"""True-point check of herm4 (swap-reduced two-letter relaxation): random orthonormal u, v; X_c = twirl of the letter
products from the multi-replica moments; the parametrisation must reproduce X_c exactly (this also checks the built-in
swap identities), every PSD block must be PSD, moments and identities exact, objective = q_n(C) computed directly.
usage: python validate_herm4.py n [cons]"""
import sys, time, numpy as np
from alg4 import build_algebras
from herm4 import Herm4, FAMS, moments_multi, ALL_CONS
from validate import zfull_from_moments, qn_direct

n = int(sys.argv[1]); cons = sys.argv[2].split(',') if len(sys.argv) > 2 else ALL_CONS
algs = build_algebras()
rel = Herm4(n, algs, cons)
rng = np.random.default_rng(23 + n)
for trial in range(2):
    D = 3 ** n
    X = rng.standard_normal((D, 2)) + (1j * rng.standard_normal((D, 2)) if trial == 0 else 0)
    Q, _ = np.linalg.qr(X); u = Q[:, :1]; v = Q[:, 1:2]
    t0 = time.time()
    zf = {c: zfull_from_moments(moments_multi([u if ch == 'u' else v for ch in c], n), n, algs[0]) for c in FAMS}
    x, res = rel.x_from_zfulls(zf)
    for (cu, cv) in [(1.0, -1.0), (1.0, 0.25), (1.0, -0.5)]:
        mins, mom, eqr, obj = rel.evaluate(x, cu, cv)
        C = cu * (u @ u.conj().T) + cv * (v @ v.conj().T)
        print(f"n={n} trial {trial} ({'complex' if trial == 0 else 'real'}) C=({cu},{cv}): reconstruction {res:.1e}, moment res {mom:.1e}, "
              f"identity res {eqr:.1e}, objective {obj:.12f} vs q_n(C) {qn_direct(C, n):.12f}; "
              f"min eig {', '.join(f'{k}:{v_:.1e}' for k, v_ in mins.items())}; {time.time() - t0:.1f}s", flush=True)
