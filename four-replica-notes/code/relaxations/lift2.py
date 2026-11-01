"""One certificate for every ratio: 2x2 lifting of the letter-resolved Hermitian relaxation (herm4).

The scalar relaxation (fixed C = cu P_u + cv P_v) is   min c.x  s.t.  A x + s = b, s in K   (Clarabel form), with
c = cu^2 m_uu + cv^2 m_vv + 2 cu cv m_uv.  Put M(x) = [[m_uu.x, m_uv.x], [m_uv.x, m_vv.x]].  The endpoint for every
Hermitian rank-two C (every ratio and sign) is M(x) >= 0 on the true points.

Lifted problem: variables X = (X11, X12, X22) in R^N (x) Sym(2) and E = [[e11, e12], [e12, e22]]:
    zero rows       A_eq X_ab = b_eq E_ab                          (ab = 11, 12, 22)
    every cone row  [[b_r E11 - A_r X11, ...]] : a scalar nonneg row becomes a 2x2 PSD block, an m x m PSD block
                    becomes the 2m x 2m block [[S(X11), S(X12)], [S(X12), S(X22)]] with S(Y) = mat(b E - A Y)
    E >= 0, Tr E = 1
    minimise  m_uu.X11 + 2 m_uv.X12 + m_vv.X22 .
For every feasible x of the scalar problem and every unit v, (X, E) = (x (x) v v^T, v v^T) is feasible with value
v^T M(x) v.  So lifted value >= 0 certifies M(x) >= 0 on the whole scalar feasible set, i.e. q_n(C) >= 0 for every
Hermitian rank-two C at once (a matrix-valued dual certificate), given the validity of the scalar relaxation.
usage: python lift2.py n [tol]"""
import sys, os, time, numpy as np, scipy.sparse as sp, clarabel
from alg4 import build_algebras
from herm4 import Herm4, ALL_CONS

def triu_index(m):
    """position of (i, j), i <= j, in Clarabel's column-major packed upper triangle"""
    return {(i, j): j * (j + 1) // 2 + i for j in range(m) for i in range(j + 1)}

def lift(A, b, cones, N):
    A = sp.csr_matrix(A)
    NX = 3 * N + 3                      # X11, X12, X22, e11, e12, e22
    off = {'11': 0, '12': N, '22': 2 * N}
    eidx = {'11': 3 * N, '12': 3 * N + 1, '22': 3 * N + 2}
    R, Cc, V, bb, cones2 = [], [], [], [], []
    r_out = 0
    def add_row(entries, rhs):
        """entries: list of (col, val); row of A2 with b2 = rhs"""
        nonlocal r_out
        for c, v in entries:
            R.append(r_out); Cc.append(c); V.append(v)
        bb.append(rhs); r_out += 1
    def row_entries(r, ab, scale):
        """-A_r . X_ab + b_r e_ab  as  b2 - A2 x  -> A2 entries are +A_r on X_ab and -b_r on e_ab"""
        st, en = A.indptr[r], A.indptr[r + 1]
        ent = [(off[ab] + A.indices[k], scale * A.data[k]) for k in range(st, en)]
        if b[r] != 0.0: ent.append((eidx[ab], -scale * b[r]))
        return ent
    r_in = 0
    # zero cones first (as in herm4.assemble), then nonneg, then PSD
    nzero = 0
    for cone in cones:
        name = type(cone).__name__
        if 'Zero' in name:
            dim = cone.dim if hasattr(cone, 'dim') else None
            if dim is None: raise RuntimeError('zero cone dim')
            for r in range(r_in, r_in + dim):
                for ab in ('11', '12', '22'):
                    add_row(row_entries(r, ab, 1.0), 0.0); nzero += 1
            r_in += dim
    # trace(E) = 1
    add_row([(eidx['11'], 1.0), (eidx['22'], 1.0)], 1.0); nzero += 1
    cones2.append(clarabel.ZeroConeT(nzero))
    # E >= 0  (2x2 PSD, packed [e11, sqrt2 e12, e22]; s = b - A x = -(-...)): s = [e11, sqrt2 e12, e22]
    add_row([(eidx['11'], -1.0)], 0.0); add_row([(eidx['12'], -np.sqrt(2.0))], 0.0); add_row([(eidx['22'], -1.0)], 0.0)
    cones2.append(clarabel.PSDTriangleConeT(2))
    r_in = 0
    for cone in cones:
        name = type(cone).__name__
        if 'Zero' in name:
            r_in += cone.dim; continue
        if 'Nonneg' in name:
            for r in range(r_in, r_in + cone.dim):
                # 2x2 block [[s(X11), s(X12)], [s(X12), s(X22)]], packed (0,0), (0,1)*sqrt2, (1,1)
                add_row(row_entries(r, '11', 1.0), 0.0)
                add_row(row_entries(r, '12', np.sqrt(2.0)), 0.0)
                add_row(row_entries(r, '22', 1.0), 0.0)
                cones2.append(clarabel.PSDTriangleConeT(2))
            r_in += cone.dim; continue
        # PSD triangle cone of size m
        m = cone.dim if hasattr(cone, 'dim') else None
        k = m * (m + 1) // 2
        tri = triu_index(m); tri2 = triu_index(2 * m)
        rows_small = {ij: r_in + p for ij, p in tri.items()}
        big = [None] * (2 * m * (2 * m + 1) // 2)
        for (p, q), pos in tri2.items():
            if q < m:              # both in first half -> X11 entry (p, q)
                ab, i, j = '11', p, q
            elif p >= m:           # both in second half -> X22 entry
                ab, i, j = '22', p - m, q - m
            else:                  # cross block -> X12 entry (p, q - m)
                ab, i, j = '12', p, q - m
            ii, jj = (i, j) if i <= j else (j, i)
            small_scale = 1.0 if ii == jj else np.sqrt(2.0)
            big_scale = 1.0 if p == q else np.sqrt(2.0)
            big[pos] = row_entries(rows_small[(ii, jj)], ab, big_scale / small_scale)
        for ent in big: add_row(ent, 0.0)
        cones2.append(clarabel.PSDTriangleConeT(2 * m))
        r_in += k
    assert r_in == A.shape[0]
    A2 = sp.csc_matrix((V, (R, Cc)), shape=(r_out, NX))
    return A2, np.array(bb), cones2, off, eidx

if __name__ == "__main__":
    n = int(sys.argv[1]); tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-9
    algs = build_algebras(); t0 = time.time()
    rel = Herm4(n, algs, ALL_CONS, verbose=False)
    A, b, cones = rel.assemble(); N = rel.N
    print(f"n={n}: scalar problem {A.shape}, N = {N}; build {time.time() - t0:.0f}s", flush=True)
    A2, b2, cones2, off, eidx = lift(A, b, cones, N)
    c2 = np.zeros(A2.shape[1])
    c2[off['11']:off['11'] + N] = rel.m_uu; c2[off['12']:off['12'] + N] = 2 * rel.m_uv; c2[off['22']:off['22'] + N] = rel.m_vv
    print(f"   lifted problem {A2.shape}, nnz {A2.nnz / 1e6:.1f}M (scalar nnz {sp.csc_matrix(A).nnz / 1e6:.1f}M)", flush=True)
    def k2(cs):
        tot = 0; mx = 0
        for c in cs:
            if type(c).__name__.startswith('PSD'):
                k = c.dim * (c.dim + 1) // 2; tot += k * k; mx = max(mx, c.dim)
        return tot, mx
    s1, m1 = k2(cones); s2, m2 = k2(cones2)
    print(f"   dense PSD scaling blocks: scalar sum k^2 = {s1 / 1e6:.1f}M (max m {m1}), lifted sum k^2 = {s2 / 1e6:.1f}M (max m {m2})", flush=True)
    if len(sys.argv) > 3 and sys.argv[3] == 'sizeonly': sys.exit(0)
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 400
    st.tol_gap_abs = tol; st.tol_gap_rel = tol; st.tol_feas = tol
    t1 = time.time()
    sol = clarabel.DefaultSolver(sp.csc_matrix((A2.shape[1], A2.shape[1])), c2, A2, b2, cones2, st).solve()
    x = np.array(sol.x); z = np.array(sol.z)
    E = np.array([[x[eidx['11']], x[eidx['12']]], [x[eidx['12']], x[eidx['22']]]])
    r = A2.T @ z + c2
    print(f"n={n} LIFTED (all ratios at once): status {sol.status}, iters {sol.iterations}, primal {sol.obj_val:+.6e}, "
          f"dual {-(b2 @ z):+.6e}, ||A2^T z + c2||_inf {np.abs(r).max():.1e}, E = {np.round(E, 4).tolist()}; "
          f"solve {time.time() - t1:.0f}s", flush=True)
    # sanity: scalar values at the two zero directions with the same machinery
    for (cu, cv) in [(1.0, 1.0), (1.0, -1.0), (1.0, 0.5)]:
        s1 = rel.solve(cu, cv)
        print(f"   scalar check C = ({cu}, {cv}): {s1.status}, value {s1.obj_val:+.6e}", flush=True)
