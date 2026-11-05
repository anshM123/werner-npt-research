"""Floating-point observation (no claim rests on it): the minimiser X* of the eliminated (C_u) block for the letter-2
data v of Conjecture M, and for the special direction v1 = (1, 1/s_t, ..., 1/s_t), is NOT of the separable form
x_rh = a_r + b_h off the diagonal once u >= 4 (residual of the best additive fit is of order 1e-4, far above
roundoff).  This is the pair structure created by the crossed term (rho_r rho_h/(2+g)) X_rh X_hr of the (C_u) form,
and it is why G_z has no closed form in finitely many moment sums of the U-weights (Section 3 of NOTES.md).
Also prints how many entries of the (C_u) block (scaled coordinates, gauge x_t = 0) change sign with the weights."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ex_sep import solve
from ex_csign import blocks

rng = np.random.default_rng(3)
ok = True
for u in (3, 4, 5, 6):
    st = float(np.exp(rng.uniform(-1, 1)))
    sU = list(np.exp(rng.uniform(-1.5, 1.5, u)))
    S = sum(sU)
    beta = 3 * (2 * st + S) / (4 * st * (st + S))
    for name, vz in (('v', np.array([1.0] + [beta] * u)), ('v1', np.array([1.0] + [1.0 / st] * u))):
        Xs, q, X = solve(st, sU, vz)
        rows, rhs = [], []
        for r in range(u):
            for h in range(u):
                if r != h:
                    a = np.zeros(2 * u); a[r] = 1; a[u + h] = 1
                    rows.append(a); rhs.append(Xs[r, h])
        A = np.array(rows)
        sol = np.linalg.lstsq(A, np.array(rhs), rcond=None)[0]
        res = np.abs(np.array(rhs) - A @ sol).max()
        asym = np.abs(Xs - Xs.T).max()
        print(f'u = {u}, {name}: additive-fit residual {res:.2e}, asymmetry of x* {asym:.2e}')
        if u >= 4:
            ok &= res > 1e-6
mixed = 0
tot = 0
SC = []
for it in range(200):
    stt = np.exp(rng.uniform(-4, 4)); sUU = np.exp(rng.uniform(-4, 4, 3))
    C, B, lab = blocks(stt, sUU)
    SC.append(np.sign(C))
SC = np.array(SC)
mixed = int(((SC.max(0) > 0) & (SC.min(0) < 0)).sum())
print(f'(C_u) block, u = 3, gauge x_t = 0, scaled: {mixed} of {SC.shape[1] ** 2} entries change sign over 200 weight vectors')
print('PASS' if ok else 'FAIL')
