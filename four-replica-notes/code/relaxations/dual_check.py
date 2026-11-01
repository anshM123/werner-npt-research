"""Solve the four-replica relaxation and inspect the dual certificate.
Clarabel form: min c.x  s.t.  A x + s = b, s in K.  Dual: max -b.z  s.t.  A^T z + c = 0, z in K* (= K).
For any z in K with residual r = A^T z + c, and any primal-feasible x:  c.x = -b.z + s.z + r.x  >= -b.z + r.x.
usage: python dual_check.py n cons [tol] [rank]"""
import sys, time, numpy as np, scipy.sparse as sp
from alg4 import build_algebras
from rep4 import Relax, unsvec

def cone_blocks(rel):
    """list of (kind, size) in row order matching rel.Amat"""
    import clarabel
    out = []
    for c in rel.cones:
        name = type(c).__name__
        if 'Zero' in name: out.append(('zero', c.dim if hasattr(c, 'dim') else None))
        elif 'Nonneg' in name: out.append(('nn', None))
        else: out.append(('psd', None))
    return out

if __name__ == "__main__":
    n = int(sys.argv[1]); cons = sys.argv[2].split(','); tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
    RANK = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
    algs = build_algebras()
    rel = Relax(n, RANK, algs, cons, verbose=False)
    sol = rel.solve(tol=tol, max_iter=1000)
    x = np.array(sol.x); z = np.array(sol.z); s = np.array(sol.s)
    A = rel.Amat; b = rel.bvec; c = rel.obj
    r = A.T @ z + c
    print(f"n={n} rank={RANK:g} cons={'+'.join(cons)} tol={tol:g}: status {sol.status}, iters {sol.iterations}")
    print(f"   primal obj c.x = {c @ x:.6e}  (q_n = {c @ x / RANK ** 2:.6e});  dual obj -b.z = {-(b @ z):.6e} "
          f"(q_n >= {-(b @ z) / RANK ** 2:.6e});  ||A^T z + c||_inf = {np.abs(r).max():.2e}")
    # walk the cones to get min eigenvalues of z blocks and of s blocks
    row = 0; minz = np.inf; mins = np.inf
    m_eq = len(rel.mom_rhs) + rel.eq_rows.shape[0]
    row = m_eq
    nn = 0
    for (o, m) in rel.wcones:
        if m == 1: nn += 1
    for (cn, mu, cc, M) in rel.blocks:
        if cc['m'] == 1: nn += 1
    if nn:
        minz = min(minz, z[row:row + nn].min()); mins = min(mins, s[row:row + nn].min()); row += nn
    sizes = [m for (o, m) in rel.wcones if m > 1] + [cc['m'] for (cn, mu, cc, M) in rel.blocks if cc['m'] > 1]
    for m in sizes:
        k = m * (m + 1) // 2
        minz = min(minz, np.linalg.eigvalsh(unsvec(z[row:row + k], m)).min())
        mins = min(mins, np.linalg.eigvalsh(unsvec(s[row:row + k], m)).min()); row += k
    assert row == A.shape[0]
    print(f"   min eig of dual blocks z: {minz:.2e};  min eig of primal slacks s: {mins:.2e};  moment multipliers {z[:5]}")
    np.savez(f"dual_n{n}_{'_'.join(cons)}_r{int(RANK)}.npz", x=x, z=z, s=s, c=c, b=b)
