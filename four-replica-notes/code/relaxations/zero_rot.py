"""One zero of q_n at t = +1 with a rotated two-letter basis: u, v orthogonal products on k sites (k odd) sharing a
random psi on the other sites, then (u, v) -> (u, v) R for a random 2x2 unitary R (C = |u><u| + |v><v| is unchanged).
Saves the parameter vector x of the herm4 parametrisation.  usage: python zero_rot.py n k seed"""
import sys, numpy as np
from alg4 import build_algebras
from herm4 import Herm4, FAMS, moments_multi
from validate import zfull_from_moments, qn_direct
n = int(sys.argv[1]); k = int(sys.argv[2]); seed = int(sys.argv[3])
rng = np.random.default_rng(seed)
algs = build_algebras(); rel = Herm4(n, algs, [], verbose=False)
def unit(D):
    x = rng.standard_normal(D) + 1j * rng.standard_normal(D); return x / np.linalg.norm(x)
u = np.ones(1, complex); v = np.ones(1, complex)
for s in range(k):
    Q, _ = np.linalg.qr(rng.standard_normal((3, 2)) + 1j * rng.standard_normal((3, 2)))
    u = np.kron(u, Q[:, 0]); v = np.kron(v, Q[:, 1])
if n - k > 0:
    psi = unit(3 ** (n - k)); u = np.kron(u, psi); v = np.kron(v, psi)
R, _ = np.linalg.qr(rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2)))
U = np.stack([u, v], axis=1) @ R
u2, v2 = U[:, :1], U[:, 1:2]
zf = {c: zfull_from_moments(moments_multi([u2 if ch == 'u' else v2 for ch in c], n), n, algs[0]) for c in FAMS}
x, res = rel.x_from_zfulls(zf)
C = u2 @ u2.conj().T + v2 @ v2.conj().T
mins, mom, eqr, obj = rel.evaluate(x, 1.0, 1.0)
print(f"n={n} k={k} seed={seed}: reconstruction {res:.1e}, q(C) direct {qn_direct(C, n):.1e}, objective {obj:.1e}, moment res {mom:.1e}, "
      f"min eig X {mins.get('X', np.nan):.1e}", flush=True)
np.save(f"zero_rot_n{n}_k{k}_s{seed}.npy", x)
