"""Search step for a rigorous SF-form certificate of a two-atom target (then run verify_bm_cert.py).
The target is X = G(A1, B1) + G(A2, B2) with the four filters saved as exact float64 arrays.  We fit
T = X/Tr X - (eps/9) I by ONE atom + products (fit_robust, variable projection to machine precision); the remainder
Z = X - sum w_i G_i is then close to (eps Tr X / 9) I and the verifier certifies it separable (Gurvits-Barnum ball).
usage: python cert_search.py source.npz index tag eps
  source.npz: advmu_*.npz (mode, params, tau) -> target normalised(G^(1) + w G^(2)), rescaled into G(A1,B1) + G(A2,B2)"""
import sys, numpy as np
from bm_core import G, RHO, ptmin, D, N
from bm_fit import fit_robust, norm_t
from ppt_margin import maximize_mu
from bm_fit import candidates

src, idx, tag, EPS = sys.argv[1], int(sys.argv[2]), sys.argv[3], float(sys.argv[4])
Z = np.load(src, allow_pickle=True)
mode = str(Z['mode']); pp = Z['params'][idx]; w = float(np.exp(Z['tau'][idx]))
if mode == 'nf':
    A1, B1, A2, B2 = np.eye(D, dtype=complex), np.eye(D, dtype=complex), pp[0].copy(), pp[1].copy()
    s1 = 1 / np.sqrt(np.trace(RHO).real)
else:
    A1, B1, A2, B2 = pp[0].copy(), pp[1].copy(), pp[2].copy(), pp[3].copy()
    s1 = 1 / np.sqrt(np.trace(G(A1, B1)).real)
s2 = np.sqrt(w / np.trace(G(A2, B2)).real)
A1 = A1 * s1; A2 = A2 * s2                       # now X = G(A1,B1) + G(A2,B2) (exact floats from here on)
np.savez(f"cert_target_{tag}.npz", A1=A1, B1=B1, A2=A2, B2=B2)
X = G(A1, B1) + G(A2, B2)
Xn = norm_t(X)
print(f"target {tag} ({src}[{idx}], mode {mode}): PT-min {ptmin(Xn):+.4f}, lambda_min {np.linalg.eigvalsh(Xn)[0]:.4f}", flush=True)
m, E, Fm, lam = maximize_mu(Xn, candidates(Xn) + [(A1, B1), (A2, B2)], iters=300)
print(f"   PPT margin mu* >= {m:+.5f}", flush=True)
T = Xn - EPS / N * np.eye(N)
scale = np.trace(T).real
out = fit_robust(T / scale, P=80, restarts=3, rounds=40, seed=11, extra=[(A1, B1), (A2, B2), (E, Fm)])
fo = out['fit']
coef, r = fo.solve(fo.th)
print(f"   fit of X/Tr X - (eps/9) I, eps = {EPS:.0e}: residual {np.linalg.norm(r) / fo.nX:.3e} (abs {np.linalg.norm(r) * scale:.3e}), "
      f"atom weight {coef[0]:.4f}, active products {int(np.sum(coef[1:] > 0))}", flush=True)
coef = coef * scale
act = np.where(coef[1:] > 0)[0]
kinds = np.array(['f'] + ['p'] * len(act))
zf, zp = np.zeros((D, D), complex), np.zeros(D, complex)
Ef = (fo.th[:9] + 1j * fo.th[9:18]).reshape(D, D); Ff = (fo.th[18:27] + 1j * fo.th[27:36]).reshape(D, D)
np.savez(f"cert_atoms_{tag}.npz", kinds=kinds,
         P=np.array([Ef] + [zf] * len(act)), Q=np.array([Ff] + [zf] * len(act)),
         x=np.array([zp] + [fo.dict.xv[i] for i in act]), y=np.array([zp] + [fo.dict.yv[i] for i in act]),
         coef=np.concatenate([[coef[0]], coef[1:][act]]), eps=EPS, residual=np.linalg.norm(r) / fo.nX)
print(f"   saved cert_target_{tag}.npz, cert_atoms_{tag}.npz", flush=True)
