"""One zero of q_n at t = +1 from the interpolating family, with a rotated two-letter basis:
u = a (x) phi, v = a' (x) chi, phi = cos th w+ + e^{i al} sin th w-, chi = sin th w+ + e^{i al} cos th w-, where
(w+, w-) are orthogonal products on two sites sharing psi on the remaining n - 3 sites; then (u, v) -> (u, v) R.
Saves the herm4 parameter vector.  usage: python zero_interp.py n seed"""
import sys, numpy as np
from alg4 import build_algebras
from herm4 import Herm4, FAMS, moments_multi
from validate import zfull_from_moments, qn_direct
n = int(sys.argv[1]); seed = int(sys.argv[2]); rng = np.random.default_rng(seed)
algs = build_algebras(); rel = Herm4(n, algs, [], verbose=False)
def unit(D):
    x = rng.standard_normal(D) + 1j * rng.standard_normal(D); return x / np.linalg.norm(x)
wp = np.ones(1, complex); wm = np.ones(1, complex)
for s in range(2):
    Qp, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2))); wp = np.kron(wp, Qp[:, 0]); wm = np.kron(wm, Qp[:, 1])
if n - 3 > 0:
    psi = unit(3 ** (n - 3)); wp = np.kron(wp, psi); wm = np.kron(wm, psi)
th = rng.uniform(0, np.pi / 2); al = rng.uniform(0, 2 * np.pi)
Qa, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
u = np.kron(Qa[:, 0], np.cos(th) * wp + np.exp(1j * al) * np.sin(th) * wm)
v = np.kron(Qa[:, 1], np.sin(th) * wp + np.exp(1j * al) * np.cos(th) * wm)
R, _ = np.linalg.qr(rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2)))
U = np.stack([u, v], axis=1) @ R
zf = {c: zfull_from_moments(moments_multi([U[:, :1] if ch == 'u' else U[:, 1:2] for ch in c], n), n, algs[0]) for c in FAMS}
x, res = rel.x_from_zfulls(zf)
C = U[:, :1] @ U[:, :1].conj().T + U[:, 1:2] @ U[:, 1:2].conj().T
mins, mom, eqr, obj = rel.evaluate(x, 1.0, 1.0)
print(f"n={n} interpolating seed={seed} theta={th:.3f}: reconstruction {res:.1e}, q(C) direct {qn_direct(C, n):.1e}, objective {obj:.1e}, "
      f"min eig X {mins.get('X', np.nan):.1e}", flush=True)
np.save(f"zero_rot_n{n}_interp_s{seed}.npy", x)
