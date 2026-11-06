import numpy as np, flint, time
from psi_exact import PsiExact, ldl_pivots
from relax2 import Psi_matrix
rng = np.random.default_rng(2)
for u in (3, 4, 6):
    wt = flint.fmpq(int(rng.integers(1, 20)), int(rng.integers(1, 20)))
    wU = [flint.fmpq(int(rng.integers(1, 30)), int(rng.integers(1, 30))) for _ in range(u)]
    t0 = time.time()
    P = PsiExact(wt, wU)
    Mx, iy = P.matrix()
    Mf = np.array([[float(Mx[i, j]) for j in range(Mx.ncols())] for i in range(Mx.nrows())])
    st = float(wt) ** 2
    sU = np.array([float(x) ** 2 for x in wU])
    M2 = Psi_matrix(st, sU, 'V2')
    k = u + 1

    def Gr(Mm):
        A = Mm[:-k, :-k]; B = Mm[:-k, -k:]; D = Mm[-k:, -k:]
        return D - B.T @ np.linalg.solve(A, B)
    g1 = Gr(Mf); g2 = Gr(M2)
    piv = ldl_pivots(Mx)
    print(u, 'max |G_rel exact - float|', np.abs(g1 - g2).max(), ' all pivots > 0:', all(p > 0 for p in piv),
          f'{time.time() - t0:.2f}s', flush=True)
