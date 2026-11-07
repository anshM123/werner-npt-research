"""Exact (integer) verification of the first-letter certificate at a = 3/4 on Sym^n (x) Sym^n.

Unnormalised symmetric basis f_alpha = sum of e_w over words w of content alpha (integer vectors), so
f_alpha = sum_i e_i (x) f_{alpha - d_i}.  For one-copy integer operators X and an operator Y on n-1 copies,
<f_a f_b | X (x) Y | f_a' f_b'> = sum_{ijkl} X[(i,j),(k,l)] <f_{a-d_i} f_{b-d_j} | Y | f_{a'-d_k} f_{b'-d_l}>.

T_m = f-matrix of (4C)^{(x)m},  4C = 3I - F - Omega Omega^T  (C = Choi matrix of P_{3/4});
U_m = f-matrix of 4^m Q_m,  Q_m = sum_{j<m} (I/4)^{(x)j} (x) (A/2) (x) C^{(x)(m-1-j)},  A = P_an = (I - F)/2,
    U_m = tr(I, U_{m-1}) + tr(I - F, T_{m-1}).
Telescoping gives Q_m + Gamma(Q_m) = C^{(x)m} - 4^{-m} I exactly (checked below as an integer identity).
Positivity: Q_n is SO(3)-invariant, so Q_n >= 0 iff its restriction to ker(J_z) is >= 0; ker J_z (in Sym (x) Sym) is
spanned by the SO(2)-invariant bihomogeneous polynomials r_u^a r_w^b s^c tau^eps u3^d w3^e (eps in {0,1}),
r = x1^2 + x2^2, s = u1w1 + u2w2, tau = u1w2 - u2w1.  The Gram matrix M = B^T U_n B is integral (B scaled),
split into parity blocks, and checked PSD by exact LDL^T over the rationals (python-flint fmpq).
"""
import sys, time
import numpy as np
from math import factorial
import flint
from bos import mindices, index_of


def one_copy_int(d=3):
    I = np.eye(d * d, dtype=np.int64)
    F = np.zeros((d * d, d * d), dtype=np.int64)
    Om = np.zeros((d * d, d * d), dtype=np.int64)
    for i in range(d):
        for j in range(d):
            F[i * d + j, j * d + i] = 1
    for i in range(d):
        for k in range(d):
            Om[i * d + i, k * d + k] = 1
    return I, F, Om


def shift_maps(m, d=3):
    """For level m: for each i, index map alpha -> index(alpha - d_i) at level m-1 (or -1)."""
    idx = index_of(m - 1, d)
    res = []
    for i in range(d):
        t = []
        for al in mindices(m, d):
            if al[i] > 0:
                be = list(al); be[i] -= 1; t.append(idx[tuple(be)])
            else:
                t.append(-1)
        res.append(np.array(t))
    return res


def tr(X, Y, m, d=3):
    """Integer transfer: level m-1 matrix Y -> level m matrix, one-copy operator X (9x9 integer)."""
    Nm1 = len(mindices(m - 1, d)); N = len(mindices(m, d))
    sh = shift_maps(m, d)
    out = np.zeros((N * N, N * N), dtype=np.int64)
    pair_idx = {}; pair_ok = {}
    for i in range(d):
        for j in range(d):
            ti, tj = sh[i], sh[j]
            ok = (ti[:, None] >= 0) & (tj[None, :] >= 0)
            pidx = (np.maximum(ti, 0)[:, None] * Nm1 + np.maximum(tj, 0)[None, :]).reshape(-1)
            pair_idx[i, j] = pidx; pair_ok[i, j] = ok.reshape(-1)
    for i in range(d):
        for j in range(d):
            for k in range(d):
                for l in range(d):
                    c = int(X[i * d + j, k * d + l])
                    if c == 0:
                        continue
                    r = pair_idx[i, j]; cidx = pair_idx[k, l]
                    blk = Y[np.ix_(r, cidx)] * (pair_ok[i, j][:, None] & pair_ok[k, l][None, :])
                    out += c * blk
    return out


def build_UT(n):
    I, F, Om = one_copy_int()
    C4 = 3 * I - F - Om
    T = np.array([[1]], dtype=np.int64); U = np.array([[0]], dtype=np.int64); E = np.array([[1]], dtype=np.int64)
    lim = 2 ** 63 // 64   # each output entry is a sum of at most 45 terms |X_ijkl| * |Y| (sum |X| <= 45)
    for m in range(1, n + 1):
        for Y in (U, T, E):
            assert int(np.abs(Y).max()) < lim, ('int64 overflow risk at level', m)
        U = tr(I, U, m) + tr(I - F, T, m)
        T = tr(C4, T, m)
        E = tr(I, E, m)
    print(f"   max|U| = {int(np.abs(U).max()):.3e}, max|T| = {int(np.abs(T).max()):.3e} (int64 bound check passed at every level)")
    return U, T, E


def check_float(n, U):
    """Overflow guard: recompute in float64 and compare (entries below 2^53 are exact in float)."""
    I, F, Om = [x.astype(float) for x in one_copy_int()]
    C4 = 3 * I - F - Om
    T = np.array([[1.0]]); Uf = np.array([[0.0]])
    for m in range(1, n + 1):
        Uf = trf(I, Uf, m) + trf(I - F, T, m)
        T = trf(C4, T, m)
    return np.abs(Uf - U.astype(float)).max(), np.abs(Uf).max()


def trf(X, Y, m, d=3):
    Nm1 = len(mindices(m - 1, d)); N = len(mindices(m, d))
    sh = shift_maps(m, d)
    out = np.zeros((N * N, N * N))
    pair_idx = {}; pair_ok = {}
    for i in range(d):
        for j in range(d):
            ti, tj = sh[i], sh[j]
            ok = (ti[:, None] >= 0) & (tj[None, :] >= 0)
            pair_idx[i, j] = (np.maximum(ti, 0)[:, None] * Nm1 + np.maximum(tj, 0)[None, :]).reshape(-1)
            pair_ok[i, j] = ok.reshape(-1)
    for i in range(d):
        for j in range(d):
            for k in range(d):
                for l in range(d):
                    c = X[i * d + j, k * d + l]
                    if c == 0:
                        continue
                    blk = Y[np.ix_(pair_idx[i, j], pair_idx[k, l])] * (pair_ok[i, j][:, None] & pair_ok[k, l][None, :])
                    out += c * blk
    return out


def poly_mul(p, q):
    out = {}
    for k1, v1 in p.items():
        for k2, v2 in q.items():
            k = tuple(a + b for a, b in zip(k1, k2))
            out[k] = out.get(k, 0) + v1 * v2
    return {k: v for k, v in out.items() if v != 0}


def poly_pow(p, e):
    out = {(0,) * 6: 1}
    for _ in range(e):
        out = poly_mul(out, p)
    return out


def invariant_basis(n):
    """SO(2)-invariant bihomogeneous polynomials of bidegree (n,n) in u=(u1,u2,u3), w=(w1,w2,w3).
    Exponent tuples (u1,u2,u3,w1,w2,w3).  Returns list of (label, polydict)."""
    ru = {(2, 0, 0, 0, 0, 0): 1, (0, 2, 0, 0, 0, 0): 1}
    rw = {(0, 0, 0, 2, 0, 0): 1, (0, 0, 0, 0, 2, 0): 1}
    s = {(1, 0, 0, 1, 0, 0): 1, (0, 1, 0, 0, 1, 0): 1}
    tau = {(1, 0, 0, 0, 1, 0): 1, (0, 1, 0, 1, 0, 0): -1}
    u3 = {(0, 0, 1, 0, 0, 0): 1}; w3 = {(0, 0, 0, 0, 0, 1): 1}
    out = []
    for eps in (0, 1):
        for c in range(n + 1):
            for a in range((n - c - eps) // 2 + 1):
                dd = n - 2 * a - c - eps
                if dd < 0:
                    continue
                for b in range((n - c - eps) // 2 + 1):
                    ee = n - 2 * b - c - eps
                    if ee < 0:
                        continue
                    p = {(0,) * 6: 1}
                    for fac, e in ((ru, a), (rw, b), (s, c), (tau, eps), (u3, dd), (w3, ee)):
                        p = poly_mul(p, poly_pow(fac, e))
                    out.append(((a, b, c, eps, dd, ee), p))
    return out


def ldl_psd(M):
    """Exact LDL^T PSD test of a symmetric fmpq_mat (list of lists of fmpq). Returns (is_psd, rank, n_zero_pivots)."""
    n = len(M)
    A = [row[:] for row in M]
    rank = 0; zeros = 0
    for k in range(n):
        p = A[k][k]
        if p < 0:
            return False, rank, zeros, k
        if p == 0:
            for j in range(k + 1, n):
                if A[k][j] != 0:
                    return False, rank, zeros, k
            zeros += 1
            continue
        rank += 1
        rowk = A[k]
        for i in range(k + 1, n):
            if rowk[i] == 0:
                continue
            f = rowk[i] / p
            Ai = A[i]
            for j in range(i, n):
                if rowk[j] != 0:
                    Ai[j] -= f * rowk[j]
            # keep symmetric storage: update lower part lazily (we only read A[i][j] for j >= i)
    return True, rank, zeros, None


def main(n):
    t0 = time.time()
    N = len(mindices(n)); idx = index_of(n)
    U, T, E = build_UT(n)
    t1 = time.time()
    # integer identity: U + Gamma(U) = T - E  (since 4^n (Q + Gamma Q) = (4C)^{(x)n} - I)
    G = lambda X: X.reshape(N, N, N, N).transpose(0, 3, 2, 1).reshape(N * N, N * N)
    ident = np.array_equal(U + G(U), T - E)
    if n <= 8:
        ferr, fmax = check_float(n, U)
    else:
        ferr, fmax = float('nan'), float(np.abs(U).max())
    print(f"n={n}: N={N}, built U,T in {t1 - t0:.1f}s; identity U + Gamma U = T - E exact: {ident}; "
          f"float cross-check max diff {ferr:.1e} (max |U| {fmax:.2e}, int64 safe: {fmax < 2**62})", flush=True)
    # invariant basis -> f-coordinates c_{ab} = p_{ab} a! b! / (n!)^2 ; scale by (n!)^2 to integers
    basis = invariant_basis(n)
    cols = []; labels = []
    for lab, p in basis:
        v = [0] * (N * N)
        for e, coef in p.items():
            al = (e[0], e[1], e[2]); be = (e[3], e[4], e[5])
            ia = idx[al]; ib = idx[be]
            fa = factorial(e[0]) * factorial(e[1]) * factorial(e[2])
            fb = factorial(e[3]) * factorial(e[4]) * factorial(e[5])
            v[ia * N + ib] += coef * fa * fb
        cols.append(v); labels.append(lab)
    d = len(cols)
    print(f"   invariant basis size {d}", flush=True)
    Bm = flint.fmpz_mat(N * N, d, [cols[j][i] for i in range(N * N) for j in range(d)])
    Um = flint.fmpz_mat(N * N, N * N, [int(x) for x in U.reshape(-1)])
    M = Bm.transpose() * (Um * Bm)
    # parity blocks: (eps, (d+e) mod 2, c+? ) -- reflection u3->-u3,w3->-w3 gives (dd+ee) parity; u2->-u2,w2->-w2 gives eps
    groups = {}
    for j, lab in enumerate(labels):
        key = (lab[3], (lab[4] + lab[5]) % 2)
        groups.setdefault(key, []).append(j)
    # check off-block entries vanish
    offmax = 0
    for j1 in range(d):
        for j2 in range(d):
            k1 = (labels[j1][3], (labels[j1][4] + labels[j1][5]) % 2)
            k2 = (labels[j2][3], (labels[j2][4] + labels[j2][5]) % 2)
            if k1 != k2 and M[j1, j2] != 0:
                offmax += 1
    print(f"   nonzero off-parity-block entries: {offmax}", flush=True)
    allpsd = True; totzero = 0; totrank = 0
    for key, js in sorted(groups.items()):
        Mb = [[flint.fmpq(M[a, b]) for b in js] for a in js]
        ok, rank, zeros, bad = ldl_psd(Mb)
        allpsd &= ok; totzero += (len(js) - rank); totrank += rank
        print(f"   block {key}: size {len(js)}, PSD {ok}, rank {rank}, kernel {len(js) - rank}", flush=True)
    print(f"   RESULT n={n}: Q_n PSD (exact) = {allpsd}; kernel dim on ker J_z = {totzero} (expected n+1 = {n + 1}); "
          f"total {time.time() - t0:.1f}s", flush=True)


if __name__ == '__main__':
    for n in [int(s) for s in sys.argv[1:]]:
        main(n)
