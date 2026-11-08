"""Solve the reduced real level-2 SDP from saved data (sym_real_l2.py ... save).
usage: python sym_solve.py n a_num a_den mode [solver] [tfix]
   mode maxt: maximise t;  center: fix t = tfix, maximise s with every block >= s I.
Reports the equality residual and the least eigenvalue of every block, and saves the solution."""
import sys
import time
import numpy as np
import cvxpy as cp
from fractions import Fraction

T0 = time.time()
n = int(sys.argv[1]); anum, aden = int(sys.argv[2]), int(sys.argv[3]); a = anum / aden
mode = sys.argv[4]; solver = sys.argv[5] if len(sys.argv) > 5 else 'CLARABEL'
tfix = float(Fraction(sys.argv[6])) if len(sys.argv) > 6 else 0.0
Dd = np.load(f'symdata_n{n}.npz')
NO = int(Dd['NO']); MS = Dd['MS']; orb = Dd['orb_of_ms']; oc = Dd['oc']; Si = Dd['Si']
mults, dims = Dd['mults'], Dd['dims']
nb = len(dims)
G = [Dd[f'G{i}'] for i in range(nb)]
# target orbit sums for C_a = a I - (F + Omega)/4 in diagram coordinates (diagram 0 = identity; indices as in sym_real_l2)
import itertools
K = 3
def matchings(points):
    if not points:
        yield []; return
    p = points[0]
    for i in range(1, len(points)):
        rest = points[1:i] + points[i + 1:]
        for mm in matchings(rest):
            yield [(p, points[i])] + mm
def canon(mm):
    return tuple(sorted(tuple(sorted(p)) for p in mm))
IDENT = canon([(s, s + K) for s in range(K)])
DIAG = [IDENT] + [mm for mm in sorted(set(canon(mm) for mm in matchings(list(range(2 * K))))) if mm != IDENT]
IDX = {mm: i for i, mm in enumerate(DIAG)}
cvec = np.zeros(len(DIAG))
cvec[IDX[canon(((0, 3), (1, 4), (2, 5)))]] += a
cvec[IDX[canon(((0, 4), (1, 3), (2, 5)))]] -= 0.25
cvec[IDX[canon(((0, 1), (3, 4), (2, 5)))]] -= 0.25
Sc = np.bincount(orb, weights=oc * np.prod(cvec[MS], axis=1), minlength=NO)
t = cp.Variable() if mode == 'maxt' else tfix
s = cp.Variable() if mode == 'center' else 0.0
Zs, expr, cons = [], 0, []
for i in range(nb):
    d = int(dims[i])
    Zb = cp.Variable((d, d), symmetric=True)
    Zs.append(Zb)
    cons.append(Zb >> s * np.eye(d))
    expr = expr + mults[i] * (G[i].reshape(NO, d * d) @ cp.reshape(Zb, (d * d,), order='C'))
cons.append(expr == Sc - t * Si)
prob = cp.Problem(cp.Maximize(t if mode == 'maxt' else s), cons)
if solver == 'SCS':
    prob.solve(solver=cp.SCS, eps=1e-10, max_iters=500000)
else:
    prob.solve(solver=cp.CLARABEL, max_iter=1000, tol_gap_abs=1e-12, tol_gap_rel=1e-12, tol_feas=1e-12)
tv = t.value if mode == 'maxt' else tfix
Zv = [Z.value for Z in Zs]
lhs = sum(mults[i] * (G[i].reshape(NO, -1) @ Zv[i].reshape(-1)) for i in range(nb))
res = np.abs(lhs - (Sc - tv * Si)).max()
mins = [np.linalg.eigvalsh(Z)[0] for Z in Zv]
print(f'n = {n}, a = {a}, mode {mode}, solver {solver}: status {prob.status}, t = {tv:.9e}, '
      f's = {(s.value if mode == "center" else 0):.3e}, max residual {res:.2e}, least block eigenvalue {min(mins):.3e}, '
      f'max |Sc| {np.abs(Sc).max():.2e}  [{time.time() - T0:.0f}s]', flush=True)
np.savez(f'symsol_n{n}_{anum}_{aden}_{mode}.npz', t=tv, **{f'Z{i}': Zv[i] for i in range(nb)})
