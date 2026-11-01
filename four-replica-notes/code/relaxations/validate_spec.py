"""True-point check of the relaxation for a rank-2 PSD operator omega with support spectrum (lam, 1 - lam):
random omega = U diag(lam, 1-lam) U^dagger, Z = twirl(omega^(x)4) from the moments, all constraint blocks must be PSD,
moments and equalities exact, objective = q_n(omega).  usage: python validate_spec.py n lam [cons]"""
import sys, time, numpy as np
from alg4 import build_algebras
from rep4 import Relax
from validate import moments, zfull_from_moments, qn_direct

n = int(sys.argv[1]); lam = float(sys.argv[2])
cons = sys.argv[3].split(',') if len(sys.argv) > 3 else ['LE', 'PT1', 'PT2', 'PT1LE', 'PT2LE', 'QTR1', 'QTR2', 'CUB']
algs = build_algebras()
rel = Relax(n, 2, algs, cons, verbose=False, spec=[lam, 1 - lam])
rng = np.random.default_rng(7 + n)
for trial in range(2):
    D = 3 ** n
    X = rng.standard_normal((D, 2)) + (1j * rng.standard_normal((D, 2)) if trial == 0 else 0)
    U, _ = np.linalg.qr(X)
    om = U @ np.diag([lam, 1 - lam]) @ U.conj().T
    t0 = time.time()
    Zf = zfull_from_moments(moments(om, n), n, algs[0])
    w, res = rel.w_from_zfull(Zf)
    mins, mres, eres, objv = rel.evaluate(w)
    q = qn_direct(om, n)
    print(f"n={n} lam={lam} trial {trial}: reconstruction {res:.1e}, moment res {mres:.1e}, equality res {eres:.1e}, "
          f"objective {objv:.12f} vs q_n(omega) {q:.12f}; min eig {', '.join(f'{k}:{v:.1e}' for k, v in mins.items())}; "
          f"{time.time() - t0:.1f}s", flush=True)
