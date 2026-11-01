"""Exact symmetry-adapted compressions of the covariance blocks.

For a covariance block (xy, mu) the Gram-weighted form B(x) on Q^{D_mu} (one copy per (r, s) entry) is invariant
under z -> R_g z for g in G = T x Stab(mu): T = {id, tau = slot swap (01)} for xy in {uu, vv} (acting as
(x)_s rho^{B1}_{mu_s}(tau)), T = {id} for xy = uv, and Stab(mu) the site permutations fixing mu (a Young subgroup).
For an irreducible type sigma of G and a rational primitive idempotent e_sigma (product of (1 +- tau)/2 and Young
symmetrisers), U_sigma = basis of range(R(e_sigma)) spans C^{m_sigma} (x) e S_sigma, and
B >= 0  <=>  U_sigma^T B U_sigma >= 0 for every sigma (for invariant B).  For any B, U^T B U >= 0 is implied by B >= 0,
so these compressions are valid constraints regardless of invariance; invariance only makes them lossless.
"""
import itertools, math, numpy as np, flint
from core4 import fm_kron, perm_matrix_fm, fm_den, fm_to_int, sorted_patterns, stabilizer

def partitions(a):
    if a == 0: yield (); return
    def rec(n, mx):
        if n == 0: yield (); return
        for k in range(min(n, mx), 0, -1):
            for rest in rec(n - k, k): yield (k,) + rest
    yield from rec(a, a)

def perm_sign(p):
    s = 1; seen = set()
    for i in range(len(p)):
        if i in seen: continue
        j = i; L = 0
        while j not in seen: seen.add(j); j = p[j]; L += 1
        if L % 2 == 0: s = -s
    return s

def young_symmetriser(shape, positions, lam_part, n):
    """R(a_T) R(b_T) for the row-reading tableau of lam_part on the given site positions (site permutations of n
    sites acting on the tensor space of the given shape)"""
    # tableau rows / columns
    rows = []; k = 0
    for r in lam_part: rows.append(positions[k:k + r]); k += r
    cols = [[rows[i][j] for i in range(len(rows)) if j < len(rows[i])] for j in range(lam_part[0])]
    def group_of(blocks):
        gens = [list(itertools.permutations(b)) for b in blocks]
        out = []
        for choice in itertools.product(*gens):
            p = list(range(n))
            for b, img in zip(blocks, choice):
                for src, dst in zip(b, img): p[src] = dst
            out.append(tuple(p))
        return out
    D = int(np.prod(shape))
    a = flint.fmpq_mat(D, D); b = flint.fmpq_mat(D, D)
    for p in group_of(rows): a = a + perm_matrix_fm(shape, p)
    for q in group_of(cols): b = b + perm_matrix_fm(shape, q) * perm_sign(q)
    return a * b

def colspace_int(E):
    """integer basis (columns, LLL-reduced) of the column space of an fmpq_mat"""
    R, rk = E.transpose().rref()
    if rk == 0: return None
    D = E.nrows()
    M = flint.fmpq_mat(rk, D, [R[i, j] for i in range(rk) for j in range(D)])
    d = fm_den(M)
    Z = flint.fmpz_mat(rk, D, [int(M[i, j] * d) for i in range(rk) for j in range(D)])
    # primitive rows, then LLL (rows are the basis vectors)
    rowsZ = []
    for i in range(rk):
        g = 0
        for j in range(D): g = math.gcd(g, int(Z[i, j]))
        rowsZ.append([int(Z[i, j]) // g for j in range(D)])
    Z = flint.fmpz_mat(rk, D, [v for r in rowsZ for v in r])
    Z = Z.lll()
    return np.array([[int(Z[i, j]) for i in range(rk)] for j in range(D)], dtype=object)    # D x rk

def cov_classes(B1, xy, mu):
    """list of (label, d_sigma, U (D_mu x m) object int array) for the covariance block (xy, mu)"""
    n = len(mu); g = tuple(B1.blocks[q]['f'] for q in mu); D = int(np.prod(g))
    tau = (1, 0, 2, 3)
    Rt = fm_kron([B1.blocks[q]['rho'][tau] for q in mu]) if xy in ('uu', 'vv') else None
    # Young subgroup blocks of mu
    blocks = []
    for val in sorted(set(mu)):
        blocks.append([s for s in range(n) if mu[s] == val])
    I = flint.fmpq_mat(D, D, [1 if i == j else 0 for i in range(D) for j in range(D)])
    tparts = [('+', (I + Rt) * flint.fmpq(1, 2)), ('-', (I - Rt) * flint.fmpq(1, 2))] if Rt is not None else [('', I)]
    out = []; total = 0
    from math import factorial
    def dim_irrep(part):
        # hook length formula
        a = sum(part); hooks = 1
        for i, r in enumerate(part):
            for j in range(r):
                arm = r - j - 1; leg = sum(1 for k in range(i + 1, len(part)) if part[k] > j)
                hooks *= arm + leg + 1
        return factorial(a) // hooks
    for tl, Et in tparts:
        for parts in itertools.product(*[list(partitions(len(b))) for b in blocks]):
            E = Et
            dsig = 1
            for b, part in zip(blocks, parts):
                if len(b) > 1:
                    E = E * young_symmetriser(g, b, part, n)
                dsig *= dim_irrep(part)
            U = colspace_int(E)
            if U is None: continue
            m = U.shape[1]; total += m * dsig
            out.append((f"{tl}{'|'.join(''.join(map(str, p)) for p in parts)}", dsig, U))
    assert total == D, f"isotypic dimensions do not add up for {xy}{mu}: {total} vs {D}"
    return out

if __name__ == "__main__":
    import time
    from core4 import site_algebras
    t0 = time.time()
    A, B1, A3 = site_algebras()
    tot_svec = 0; mx = 0; nb = 0; tot_raw = 0
    for xy in ['uu', 'uv', 'vv']:
        for mu in sorted_patterns(len(B1.blocks)):
            cl = cov_classes(B1, xy, mu)
            D = int(np.prod([B1.blocks[q]['f'] for q in mu]))
            tot_raw += (2 * D) * (2 * D + 1) // 2
            for lab, dsig, U in cl:
                k = 2 * U.shape[1]; tot_svec += k * (k + 1) // 2; mx = max(mx, k); nb += 1
            if xy == 'uu' and mu in [(1, 1, 1, 1), (2, 3, 3, 3), (1, 2, 3, 3)]:
                print(xy, mu, D, [(lab, dsig, U.shape[1], max(abs(int(v)) for v in U.ravel())) for lab, dsig, U in cl])
    print(f"COV: {nb} class blocks, max size {mx}, svec rows {tot_svec} (raw {tot_raw}); {time.time() - t0:.1f}s")
