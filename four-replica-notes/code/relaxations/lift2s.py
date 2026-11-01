"""Memory-lean version of lift2.py (same lifted problem), solved with SCS.

Lifted problem (see lift2.py): variables X = (X11, X12, X22) in R^N (x) Sym(2) and E = [[e11, e12], [e12, e22]] >= 0,
Tr E = 1; every cone row of the scalar relaxation A x + s = b, s in K is lifted to a 2x2 block (scalar rows) or a
2m x 2m block (m x m PSD blocks); objective m_uu.X11 + 2 m_uv.X12 + m_vv.X22.  Value >= 0 certifies
[[m_uu.x, m_uv.x], [m_uv.x, m_vv.x]] >= 0 on the whole relaxation, i.e. every Hermitian rank-two C at once.
The lifted matrix is assembled with sparse row selections (no Python lists of entries).  SCS uses the packed lower
triangle, column-major, with off-diagonals scaled by sqrt(2).
usage: python lift2s.py n [eps] [max_iters] [indirect 0/1] [solver scs|clarabel]"""
import sys, os, time, numpy as np, scipy.sparse as sp
from alg4 import build_algebras
from herm4 import Herm4, ALL_CONS

SQ2 = np.sqrt(2.0)

def clarabel_pos(m):
    """(i, j) with i <= j -> Clarabel packed position (column-major upper triangle)"""
    return {(i, j): j * (j + 1) // 2 + i for j in range(m) for i in range(j + 1)}

def scs_order(m):
    """list of (i, j), i <= j, in SCS packed order (column-major lower triangle: entry (r, c), r >= c)"""
    out = []
    for c in range(m):
        for r in range(c, m):
            out.append((c, r))          # store as (i, j) = (c, r) with i <= j
    return out

def lift(A, b, cones, N):
    """returns A2 (csc), b2, cone dict for SCS {'z': ., 's': [...]}, with column layout X11 | X12 | X22 | e11 e12 e22"""
    A = sp.csr_matrix(A)
    NX = 3 * N + 3
    E11, E12, E22 = 3 * N, 3 * N + 1, 3 * N + 2
    def place(rows, ab, scales):
        """sparse rows (len(rows) x NX): scales[k] * A[rows[k], :] placed on block ab, plus -scales[k]*b[rows[k]] on e_ab"""
        Asub = A[rows, :]
        Asub = sp.diags(scales) @ Asub
        off = {'11': 0, '12': N, '22': 2 * N}[ab]
        ecol = {'11': E11, '12': E12, '22': E22}[ab]
        blocks = [sp.csr_matrix((len(rows), off)), Asub, sp.csr_matrix((len(rows), NX - off - N))]
        M = sp.hstack(blocks, format='csr')
        bb = np.asarray(b)[rows]
        if np.any(bb != 0):
            M = M + sp.csr_matrix((-scales * bb, (np.arange(len(rows)), np.full(len(rows), ecol))), shape=(len(rows), NX))
        return M
    out_rows = []; b2 = []; zdim = 0; sdims = []
    r_in = 0
    zero_parts = []
    for cone in cones:
        if 'Zero' in type(cone).__name__:
            rows = np.arange(r_in, r_in + cone.dim)
            for ab in ('11', '12', '22'):
                zero_parts.append(place(rows, ab, np.ones(len(rows))))
            zdim += 3 * cone.dim
        r_in += cone.dim if 'Zero' in type(cone).__name__ or 'Nonneg' in type(cone).__name__ else cone.dim * (cone.dim + 1) // 2
    # trace E = 1  (as b - A x = 0 with A row = [.. 1 at e11, 1 at e22], b = 1)
    tr = sp.csr_matrix((np.array([1.0, 1.0]), (np.array([0, 0]), np.array([E11, E22]))), shape=(1, NX))
    zero_parts.append(tr); zdim += 1
    out_rows.extend(zero_parts); b2 = [np.zeros(zdim - 1), np.array([1.0])]
    # E >= 0: s = [e11, sqrt2 e12, e22] in SCS order for m = 2: (0,0), (1,0)->(0,1), (1,1): same as Clarabel for m = 2
    out_rows.append(sp.csr_matrix((np.array([-1.0, -SQ2, -1.0]), (np.arange(3), np.array([E11, E12, E22]))), shape=(3, NX)))
    b2.append(np.zeros(3)); sdims.append(2)
    r_in = 0
    for cone in cones:
        name = type(cone).__name__
        if 'Zero' in name:
            r_in += cone.dim; continue
        if 'Nonneg' in name:
            for r in range(r_in, r_in + cone.dim):
                rows = np.array([r, r, r])
                M = sp.vstack([place(rows[:1], '11', np.array([1.0])),
                               place(rows[:1], '12', np.array([SQ2])),
                               place(rows[:1], '22', np.array([1.0]))], format='csr')
                out_rows.append(M); b2.append(np.zeros(3)); sdims.append(2)
            r_in += cone.dim; continue
        m = cone.dim; k = m * (m + 1) // 2
        cpos = clarabel_pos(m)
        order2 = scs_order(2 * m)
        sel = {'11': ([], []), '12': ([], []), '22': ([], [])}
        posmap = []
        for t, (p, q) in enumerate(order2):          # p <= q in the big matrix
            if q < m:
                ab, i, j = '11', p, q
            elif p >= m:
                ab, i, j = '22', p - m, q - m
            else:
                ab, i, j = '12', p, q - m
            ii, jj = (i, j) if i <= j else (j, i)
            small_scale = 1.0 if ii == jj else SQ2
            big_scale = 1.0 if p == q else SQ2
            sel[ab][0].append(r_in + cpos[(ii, jj)]); sel[ab][1].append(big_scale / small_scale)
            posmap.append((ab, len(sel[ab][0]) - 1))
        parts = {ab: place(np.array(sel[ab][0]), ab, np.array(sel[ab][1])) for ab in sel if len(sel[ab][0])}
        # reorder rows to the SCS packed order
        stacked = sp.vstack([parts['11'], parts['12'], parts['22']], format='csr')
        base = {'11': 0, '12': len(sel['11'][0]), '22': len(sel['11'][0]) + len(sel['12'][0])}
        perm = np.array([base[ab] + idx for ab, idx in posmap])
        out_rows.append(stacked[perm, :]); b2.append(np.zeros(len(perm))); sdims.append(2 * m)
        r_in += k
    assert r_in == A.shape[0]
    A2 = sp.vstack(out_rows, format='csc')
    return A2, np.concatenate(b2), {'z': zdim, 's': sdims}, (E11, E12, E22)

if __name__ == "__main__":
    n = int(sys.argv[1]); eps = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-6
    max_iters = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
    indirect = bool(int(sys.argv[4])) if len(sys.argv) > 4 else True
    solver = sys.argv[5] if len(sys.argv) > 5 else 'scs'
    algs = build_algebras(); t0 = time.time()
    rel = Herm4(n, algs, ALL_CONS, verbose=False)
    A, b, cones = rel.assemble(); N = rel.N
    print(f"n={n}: scalar {A.shape}, N = {N}; build {time.time() - t0:.0f}s", flush=True)
    t1 = time.time()
    A2, b2, cone2, (E11, E12, E22) = lift(A, b, cones, N)
    del A
    c2 = np.zeros(A2.shape[1]); c2[:N] = rel.m_uu; c2[N:2 * N] = 2 * rel.m_uv; c2[2 * N:3 * N] = rel.m_vv
    print(f"   lifted {A2.shape}, nnz {A2.nnz / 1e6:.1f}M, PSD blocks {len(cone2['s'])} (max {max(cone2['s'])}); "
          f"lift {time.time() - t1:.0f}s", flush=True)
    t2 = time.time()
    if solver == 'scs':
        import scs
        data = dict(A=A2, b=b2, c=c2)
        ls = scs.LinearSolver.CPU_INDIRECT if indirect else scs.LinearSolver.QDLDL
        sol = scs.SCS(data, cone2, eps_abs=eps, eps_rel=eps, max_iters=max_iters, linear_solver=ls,
                      verbose=False, acceleration_lookback=10).solve()
        x = sol['x']; y = sol['y']; info = sol['info']
        pobj = c2 @ x; dobj = -b2 @ y
        r = A2.T @ y + c2
        print(f"n={n} LIFTED via SCS (indirect={indirect}): status {info['status']}, iters {info['iter']}, primal {pobj:+.6e}, "
              f"dual {dobj:+.6e}, ||A2^T y + c||_inf {np.abs(r).max():.1e}, E = "
              f"{np.round([[x[E11], x[E12]], [x[E12], x[E22]]], 4).tolist()}; solve {time.time() - t2:.0f}s", flush=True)
    else:
        import clarabel
        raise SystemExit('use lift2.py for clarabel')
