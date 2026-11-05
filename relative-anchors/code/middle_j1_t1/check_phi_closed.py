"""Lemma P (NOTES.md, Section 1): the P-first reduction of the |T| = 1 block of (R_inf)_1 and the explicit form of Phi,
checked exactly against the Lemma B builder j1lib.py at rational weight vectors with |w| = 1, k = 3, ..., 8, every
position of t used.

Units 2^{k-1}, unrescaled, |s| = 1, sigma = x_t, U = [k] minus t, N = diag(N_r) (N_t = 1, N_r = 2 otherwise).  With the
states ordered P, then (R, Q), the block is M = [[Pi, B], [B^T, C]], and
  (a) Pi N^{-1} B = B and w^T B = 0, hence K^T M K = diag(Pi, Phi) with K = [[I, -N^{-1} B], [0, I]], Phi = C - B^T N^{-1} B;
  (b) Phi = [[D_R - (1/4) Lam P0 N^{-1} P0 Lam,  -(I + N1/2)^{-1} H + (1/4) taup (2 - sigma) Lam (E - w w_U^T) D'],
             [ . ,  2 (I + N'/2)^{-1} - taup^2 (2 - sigma) D' ((2 - sigma)/2 I - w_U w_U^T) D']],
      with P0 = I - w w^T, Lam = diag(tau_c/(1 - x_c)) (tau_c = tau on U, tau - 2 at t), D' = diag_{c in U} 1/(2 - sigma - 2x_c),
      E the embedding R^U -> R^[k], and tau, taup, N1, N', H as in Lemma T1.3 (tlib.closed);
  (c) Phi is positive definite at the point (exact LDL), so the block is >= 0 with kernel exactly p (x) |2@t> there.
Usage: python check_phi_closed.py"""
import random
from fractions import Fraction as Fr
from j1lib import J1, block_states
from tlib import closed, mul, T, inv

random.seed(23)


def unit_weights(k):
    v = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(k - 1)]
    n2 = sum(a * a for a in v)
    w = [abs(2 * a / (1 + n2)) for a in v] + [abs((1 - n2) / (1 + n2))]
    if any(a == 0 for a in w):
        return unit_weights(k)
    return w


def ldl_pd(A):
    A = [row[:] for row in A]
    n = len(A)
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True


ok = True
cnt = 0
for k in range(3, 9):
    for trial in range(3):
        w = unit_weights(k)
        assert sum(a * a for a in w) == 1
        t = (trial + k) % k
        x = [a * a for a in w]
        sig = x[t]
        st = block_states(k, (t,))
        M = J1(w).matrix(st)
        R = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 2]
        P = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
        Q = [i for i, (a, y) in enumerate(st) if a is None]
        E = R + Q
        N = [Fr(1) if r == t else Fr(2) for r in range(k)]
        Bm = [[M[P[r]][e] for e in E] for r in range(k)]
        Pi = [[M[P[r]][P[r2]] for r2 in range(k)] for r in range(k)]
        oka = all(sum(Pi[r][r2] * Bm[r2][e] / N[r2] for r2 in range(k)) == Bm[r][e] for r in range(k) for e in range(len(E)))
        oka &= all(sum(w[r] * Bm[r][e] for r in range(k)) == 0 for e in range(len(E)))
        Phi = [[M[E[a]][E[b]] - sum(Bm[r][a] * Bm[r][b] / N[r] for r in range(k)) for b in range(len(E))] for a in range(len(E))]
        # (b) explicit form
        Cf = closed(w, t)
        U = Cf['U']
        u = len(U)
        tau, taup = Cf['tau'], Cf['taup']
        Lam = [(tau - 2 if c == t else tau) / (1 - x[c]) for c in range(k)]
        Dp = [1 / (2 - sig - 2 * x[c]) for c in U]
        P0 = [[(1 if i == j else 0) - w[i] * w[j] for j in range(k)] for i in range(k)]
        Ninv = [[(1 / N[i] if i == j else 0) for j in range(k)] for i in range(k)]
        PNP = mul(mul(P0, Ninv), P0)
        A1i = inv(Cf['A1'], Fr(1))
        A1iH = mul(A1i, Cf['H'])
        A2i = inv(Cf['A2'], Fr(1))
        RR = [[(N[i] if i == j else 0) - Lam[i] * PNP[i][j] * Lam[j] / 4 for j in range(k)] for i in range(k)]
        RQ = [[-A1iH[i][j] + taup * (2 - sig) * Lam[i] * ((1 if i == U[j] else 0) - w[i] * w[U[j]]) * Dp[j] / 4
               for j in range(u)] for i in range(k)]
        QQ = [[2 * A2i[i][j] - taup ** 2 * (2 - sig) * Dp[i] * ((2 - sig) / 2 * (1 if i == j else 0) - w[U[i]] * w[U[j]]) * Dp[j]
               for j in range(u)] for i in range(u)]
        Phic = [RR[i] + RQ[i] for i in range(k)] + [T(RQ)[i] + QQ[i] for i in range(u)]
        okb = Phic == Phi
        okc = ldl_pd(Phi)
        cnt += 1
        if not (oka and okb and okc):
            print(f'FAIL k={k} t={t}: (a) {oka}, (b) {okb}, (c) {okc}', flush=True)
        ok &= oka and okb and okc
    print(f'k = {k}: P-first reduction, explicit Phi and Phi > 0, exact: {ok}', flush=True)
print(f'{cnt} blocks checked')
print('PASS' if ok else 'FAIL')
