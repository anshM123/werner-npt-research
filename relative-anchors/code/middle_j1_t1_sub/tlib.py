"""Closed form of the type I block T = {t} of (R_inf)_1 (size 3k - 1), after Lemma T1.3 of middle_j1/NOTES.md,
generic in the number type (Fractions, floats or numpy floats), and the P-first reduction Phi.

Units 2^{k-1}, unrescaled, |s| = 1 (x = s), sigma = x_t, U = [k] minus t, N_r = 2 (r in U), N_t = 1.
States in the order R_c (c in [k]), P_r (r in [k]), Q_c (c in U).

P-first reduction.  The block on the P-states is Pi = diag(N) - (N w)(N w)^T/(2 - sigma), which is >= 0 with kernel
span(w), and the coupling B = [Y_RP^T, Y_PQ] of P to (R, Q) has columns orthogonal to w.  For b orthogonal to w,
min over pi of pi^T Pi pi - 2 pi^T b equals -b^T diag(1/N) b.  Hence the block is >= 0 with kernel exactly
p (x) |2@t> = (0, w, 0) if and only if

    Phi = C - B^T diag(1/N) B  is positive definite,  C = the block on the R- and Q-states.

closed(w, t) returns the blocks; phi(w, t) returns Phi on (R, Q) (size 2k - 1)."""
from fractions import Fraction as Fr


def mat(n, m, f):
    return [[f(i, j) for j in range(m)] for i in range(n)]


def mul(A, B):
    return [[sum(A[i][l] * B[l][j] for l in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def T(A):
    return [list(r) for r in zip(*A)]


def addm(A, B, c=1):
    return [[a + c * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def inv(M, one):
    n = len(M)
    A = [list(row) + [one if i == j else 0 * one for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        r = max(range(c, n), key=lambda i: abs(A[i][c]))
        A[c], A[r] = A[r], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * yv for x, yv in zip(A[i], A[c])]
    return [row[n:] for row in A]


def closed(w, t, one=Fr(1)):
    """w with sum of squares 1 (the caller ensures it)."""
    k = len(w)
    x = [v * v for v in w]
    U = [c for c in range(k) if c != t]
    sig = x[t]
    N = [one if r == t else 2 * one for r in range(k)]
    g = [xc / (2 * (1 - xc)) for xc in x]
    gam = sum(g)
    tau = (1 + g[t]) / (1 + gam / 2)
    gamp = sum(x[u] / (2 - sig - 2 * x[u]) for u in U)
    taup = 1 / (1 + gamp / 2)
    W = w
    zero = 0 * one
    DR = mat(k, k, lambda i, j: N[i] if i == j else zero)
    taus = [tau - 2 if c == t else tau for c in range(k)]
    YRP = mat(k, k, lambda c, r: taus[c] / (2 * (1 - x[c])) * ((1 if r == c else 0) - W[r] * W[c]))
    Pi = mat(k, k, lambda r, rp: (N[r] if r == rp else 0) - N[r] * N[rp] * W[r] * W[rp] / (2 - sig))
    YPQ = mat(k, len(U), lambda r, j: taup / (2 - sig - 2 * x[U[j]]) * (N[r] * W[r] * W[U[j]] - (2 - sig) * (1 if r == U[j] else 0)))
    N1 = mat(k, k, lambda c, cp: sum(x[u] / (2 * (1 - x[c] - x[u])) for u in range(k) if u != c) if c == cp
             else W[c] * W[cp] / (2 * (1 - x[c] - x[cp])))
    H = mat(k, len(U), lambda c, j: ((1 + sig / (2 * (1 - sig - x[U[j]]))) if c == U[j] else zero) if c != t
            else W[t] * W[U[j]] / (2 * (1 - sig - x[U[j]])))
    I_k = mat(k, k, lambda i, j: one if i == j else zero)
    I_u = mat(len(U), len(U), lambda i, j: one if i == j else zero)
    A1 = addm(I_k, N1, Fr(1, 2))
    Np = mat(len(U), len(U), lambda i, j: sum(x[u] / (2 - sig - 2 * x[U[i]] - 2 * x[u]) for u in U if u != U[i]) if i == j
             else W[U[i]] * W[U[j]] / (2 - sig - 2 * x[U[i]] - 2 * x[U[j]]))
    A2 = addm(I_u, Np, Fr(1, 2))
    YRQ = [[-v for v in row] for row in mul(inv(A1, one), H)]
    YQQ = [[2 * v for v in row] for row in inv(A2, one)]
    return dict(DR=DR, YRP=YRP, Pi=Pi, YPQ=YPQ, YRQ=YRQ, YQQ=YQQ, A1=A1, A2=A2, H=H, U=U, N=N, N1=N1, Np=Np,
                tau=tau, taup=taup)


def full_block(Cf):
    """the whole block in the order R, P, Q from the closed form."""
    DR, YRP, Pi, YPQ, YRQ, YQQ = (Cf[n] for n in ('DR', 'YRP', 'Pi', 'YPQ', 'YRQ', 'YQQ'))
    rows = []
    for i in range(len(DR)):
        rows.append(DR[i] + YRP[i] + YRQ[i])
    YRPt, YPQt, YRQt = T(YRP), T(YPQ), T(YRQ)
    for i in range(len(Pi)):
        rows.append(YRPt[i] + Pi[i] + YPQ[i])
    for i in range(len(YQQ)):
        rows.append(YRQt[i] + YPQt[i] + YQQ[i])
    return rows


def phi(w, t, one=Fr(1)):
    """Phi = C - B^T diag(1/N) B on the states (R, Q), from the closed form."""
    Cf = closed(w, t, one)
    k = len(w)
    N = Cf['N']
    Bt = [Cf['YRP'][c] for c in range(k)] + T(Cf['YPQ'])          # rows: R then Q; columns: P
    C = [Cf['DR'][i] + Cf['YRQ'][i] for i in range(k)] + [T(Cf['YRQ'])[i] + Cf['YQQ'][i] for i in range(len(Cf['U']))]
    n = len(C)
    return [[C[i][j] - sum(Bt[i][r] * Bt[j][r] / N[r] for r in range(k)) for j in range(n)] for i in range(n)]
