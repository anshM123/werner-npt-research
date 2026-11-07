"""Exact verification of the GREEDY certificate at a = 3/4 on Sym^n (x) Sym^n.

C = C_{3/4} = I/4 + A/2 + B/2,  A = P_an = (I - F)/2,  B = A^Gamma = (I - Omega Omega^T)/2.
For a word w in {I, A, B}^n with k I's, s A's, t B's, weight c(w) = 4^{-k} 2^{-(s+t)}.  Greedy split:
    Q_n = (Pi x Pi) [ sum_{w: s > t} c(w) w  +  (1/2) sum_{w: s = t >= 1} c(w) w ] (Pi x Pi).
Since Gamma(A) = B, Gamma(I) = I:  Q_n + Gamma(Q_n) = K_n - 4^{-n} I  (all n; algebraic identity).
Scaled integer form: 2 * 4^n * Q_n = sum over words of (I)^k (I-F)^s (I-Omega)^t with weight 2 if s > t, 1 if s = t >= 1.
Z_m(d) := f-matrix of the sum over words of length m with s - t = d (letters I, I-F, I-Om); Z_m(-d) = Gamma Z_m(d).
Final level is done through the shifted-basis trick:  B^T tr(X, Y) B = sum_{ij} B_ij^T Y W_ij,  W_ij = sum_kl X_{ij,kl} B_kl.
"""
import sys, time
import numpy as np
import flint
from math import factorial
from bos import mindices, index_of
from exact_t2 import one_copy_int, tr, invariant_basis, ldl_psd


def gamma_int(Y, N):
    return Y.reshape(N, N, N, N).transpose(0, 3, 2, 1).reshape(N * N, N * N)


def build_Z(m0):
    """Z_m0(d) for d = 0..m0 (dict), plus E_m0 (identity f-Gram)."""
    I, F, Om = one_copy_int()
    XA = I - F; XB = I - Om
    Z = {0: np.array([[1]], dtype=np.int64)}
    E = np.array([[1]], dtype=np.int64)
    lim = 2 ** 63 // 64
    for m in range(1, m0 + 1):
        Nm1 = len(mindices(m - 1))
        get = lambda d: Z[d] if d >= 0 else (gamma_int(Z[-d], Nm1) if -d in Z else None)
        newZ = {}
        for d in range(0, m + 1):
            acc = None
            for (X, dd) in ((I, d), (XA, d - 1), (XB, d + 1)):
                Y = get(dd) if abs(dd) <= m - 1 else None
                if Y is None:
                    continue
                assert int(np.abs(Y).max()) < lim
                term = tr(X, Y, m)
                acc = term if acc is None else acc + term
            newZ[d] = acc
        E = tr(I, E, m)
        Z = newZ
    return Z, E


def basis_cols(n):
    N = len(mindices(n)); idx = index_of(n)
    basis = invariant_basis(n)
    cols = []; labels = []
    for lab, p in basis:
        v = {}
        for e, coef in p.items():
            ia = idx[(e[0], e[1], e[2])]; ib = idx[(e[3], e[4], e[5])]
            fa = factorial(e[0]) * factorial(e[1]) * factorial(e[2]); fb = factorial(e[3]) * factorial(e[4]) * factorial(e[5])
            v[(ia, ib)] = v.get((ia, ib), 0) + coef * fa * fb
        cols.append(v); labels.append(lab)
    return cols, labels


def shifted_B(n, cols, i, j):
    """B_ij as fmpz_mat on level n-1 index pairs: B_ij[(a',b'), r] = B[(a'+d_i, b'+d_j), r]."""
    Np = len(mindices(n - 1)); N = len(mindices(n))
    idxn = index_of(n)
    alp = mindices(n - 1)
    up_i = [idxn[tuple(a[t] + (1 if t == i else 0) for t in range(3))] for a in alp]
    up_j = [idxn[tuple(a[t] + (1 if t == j else 0) for t in range(3))] for a in alp]
    d = len(cols)
    rows = []
    for ap in range(Np):
        for bp in range(Np):
            key = (up_i[ap], up_j[bp])
            rows.append([cols[r].get(key, 0) for r in range(d)])
    return flint.fmpz_mat(rows)


def to_fmpz(Y):
    return flint.fmpz_mat([[int(x) for x in row] for row in Y])


def main(n):
    t0 = time.time()
    m0 = n - 1
    Z, E = build_Z(m0)
    Np = len(mindices(m0))
    G = lambda Y: gamma_int(Y, Np)
    # pieces at level m0
    P1 = sum(Z[d] for d in range(1, m0 + 1)) if m0 >= 1 else np.zeros_like(Z[0])
    P2 = sum(Z[d] for d in range(2, m0 + 1)) if m0 >= 2 else np.zeros_like(Z[0])
    Z1 = Z[1] if m0 >= 1 else np.zeros_like(Z[0])
    # 2*4^n*Q_n = tr(I, 2*P1 + Z0 - E) + tr(I-F, 2*(P1 + Z0) + Gamma Z1) + tr(I-Om, 2*P2 + Z1)
    #   (words with s>t get weight 2, s=t>=1 weight 1; the all-I word is removed via -E)
    Y_I = 2 * P1 + Z[0] - E
    Y_A = 2 * (P1 + Z[0]) + G(Z1)
    Y_B = 2 * P2 + Z1
    del Z
    t1 = time.time()
    print(f"n={n}: level {m0} built in {t1 - t0:.1f}s; max|Y| = {max(int(np.abs(Y).max()) for Y in (Y_I, Y_A, Y_B)):.3e}", flush=True)
    cols, labels = basis_cols(n)
    dB = len(cols)
    Bs = {(i, j): shifted_B(n, cols, i, j) for i in range(3) for j in range(3)}
    zero = flint.fmpz_mat(Np * Np, dB)
    Bdiag = Bs[0, 0] + Bs[1, 1] + Bs[2, 2]
    M = flint.fmpz_mat(dB, dB)
    for (Yn, kind) in ((Y_I, 'I'), (Y_A, 'A'), (Y_B, 'B')):
        Yf = to_fmpz(Yn)
        for i in range(3):
            for j in range(3):
                if kind == 'I':
                    W = Bs[i, j]
                elif kind == 'A':
                    W = Bs[i, j] - Bs[j, i]
                else:
                    W = Bs[i, j] - (Bdiag if i == j else zero)
                M += Bs[i, j].transpose() * (Yf * W)
        print(f"   accumulated {kind} part, {time.time() - t0:.1f}s", flush=True)
    sym = all(M[a, b] == M[b, a] for a in range(dB) for b in range(a + 1, dB))
    groups = {}
    for jx, lab in enumerate(labels):
        groups.setdefault((lab[3], (lab[4] + lab[5]) % 2), []).append(jx)
    off = sum(1 for a in range(dB) for b in range(dB)
              if (labels[a][3], (labels[a][4] + labels[a][5]) % 2) != (labels[b][3], (labels[b][4] + labels[b][5]) % 2) and M[a, b] != 0)
    allpsd = True; ker = 0
    for key, js in sorted(groups.items()):
        Mb = [[flint.fmpq(M[a, b]) for b in js] for a in js]
        ok, rank, zeros, bad = ldl_psd(Mb)
        allpsd &= ok; ker += len(js) - rank
        print(f"   block {key}: size {len(js)}, PSD {ok}, rank {rank}", flush=True)
    print(f"   RESULT n={n} (greedy): symmetric {sym}, off-block nonzeros {off}, Q_n PSD (exact) = {allpsd}, "
          f"kernel on ker J_z = {ker} (expected n+1 = {n + 1}); {time.time() - t0:.1f}s", flush=True)


if __name__ == '__main__':
    for n in [int(s) for s in sys.argv[1:]]:
        main(n)
