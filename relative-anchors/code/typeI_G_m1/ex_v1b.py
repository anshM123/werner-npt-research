import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ex_sep import solve
np.set_printoptions(linewidth=200, precision=12, suppress=True)
rng = np.random.default_rng(3)
for u in (3, 4, 5, 6):
    st = float(np.exp(rng.uniform(-1, 1)))
    sU = np.exp(rng.uniform(-1.5, 1.5, u))
    S = sU.sum()
    v1 = np.array([1.0] + [1.0 / st] * u)
    Xs, q, X = solve(st, list(sU), v1)
    rows, rhs = [], []
    for r in range(u):
        for h in range(u):
            if r != h:
                a = np.zeros(2 * u); a[r] = 1; a[u + h] = 1
                rows.append(a); rhs.append(Xs[r, h])
    A = np.array(rows)
    sol = np.linalg.lstsq(A, np.array(rhs), rcond=None)[0]
    res = np.abs(np.array(rhs) - A @ sol).max()
    # multiplicative: off-diagonal x_rh = f_r f_h ?
    print(f'u={u}: asym {np.abs(Xs - Xs.T).max():.2e}, additive residual {res:.2e}')
    # does x_rh (r != h) depend on a function of s_r times s_h? print x_rh vs (s_r, s_h)
    if u == 4:
        for r in range(u):
            for h in range(r + 1, u):
                print(f'   s_r={sU[r]:.4f} s_h={sU[h]:.4f} x_rh={Xs[r,h]:.12f}  a_r+a_h={sol[r]+sol[u+h]:.12f}')
        print('   diag x_hh', np.diag(Xs))
