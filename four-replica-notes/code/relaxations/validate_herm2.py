"""True-point check of herm2: random orthonormal u, v; X_c = twirl of the letter products; every block PSD,
moments and equalities exact, objective = q_n(C) for C = mu|u><u| + s nu|v><v|.  usage: python validate_herm2.py n"""
import sys, time, numpy as np
from alg4 import build_algebras
from herm2 import Herm2, FAMS, moments_multi
from validate import zfull_from_moments, qn_direct

n = int(sys.argv[1])
algs = build_algebras()
rel = Herm2(n, algs, ['SWAP', 'ORTH', 'MARG', 'PT', 'COV', 'LE2'])
rng = np.random.default_rng(11 + n)
for trial in range(2):
    D = 3 ** n
    X = rng.standard_normal((D, 2)) + (1j * rng.standard_normal((D, 2)) if trial == 0 else 0)
    Q, _ = np.linalg.qr(X); u = Q[:, :1]; v = Q[:, 1:2]
    t0 = time.time()
    zf = {}
    for c in FAMS:
        Us = [u if ch == 'u' else v for ch in c]
        zf[c] = zfull_from_moments(moments_multi(Us, n), n, algs[0])
    x, res = rel.x_from_zfulls(zf)
    for (mu, nu, s) in [(0.8, 0.6, -1.0), (0.6, 0.4, 1.0)]:
        mins, mom, eqr, obj = rel.evaluate(x, mu, nu, s)
        C = mu * (u @ u.conj().T) + s * nu * (v @ v.conj().T)
        print(f"n={n} trial {trial} (mu,nu,s)=({mu},{nu},{s:+.0f}): reconstruction {res:.1e}, moment res {mom:.1e}, eq res {eqr:.1e}, "
              f"objective {obj:.12f} vs q_n(C) {qn_direct(C, n):.12f}; min eig {', '.join(f'{k}:{v_:.1e}' for k, v_ in mins.items())}; "
              f"{time.time() - t0:.1f}s", flush=True)
