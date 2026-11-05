"""Lemma T1.3 (NOTES.md, Section 3): closed form of the type I block T = {t} of (R_inf)_1 and its inverse-free
reduction, checked exactly against the Lemma B builder j1lib at rational weights, k = 3..8.

Units 2^{k-1}, unrescaled, |s| = 1 (x = s/|s|), sigma = x_t, U = [k] minus t, N_r = 2 (r in U), N_t = 1.
  D_R = diag(N)                                   (states R_c = (2@t, 1@c))
  Y_RP[c, r] = tau_c/(2(1 - x_c)) (delta_rc - w_r w_c),   tau_c = tau (c in U), tau - 2 (c = t),
               tau = (1 + g_t)/(1 + gamma/2), g_t = sigma/(2(1 - sigma)), gamma = sum_c x_c/(2(1 - x_c))
  Pi[r, r'] = N_r delta - N_r N_r' w_r w_r'/(2 - sigma)            (states P_r = (1@r, 2@t))
  Y_PQ[r, c] = tau'/(2 - sigma - 2x_c) (N_r w_r w_c - (2 - sigma) delta_rc),  tau' = 1/(1 + gamma'/2),
               gamma' = sum_{u in U} x_u/(2 - sigma - 2x_u)          (states Q_c = (vac, 2@t 1@c), c in U)
  Y_RQ = -(I + N1/2)^{-1} H,  N1 = J^- C2^{-1} J^+ on the singles of [k]:
               N1[c,c] = sum_{u != c} x_u/(2(1 - x_c - x_u)),  N1[c,c'] = w_c w_c'/(2(1 - x_c - x_c')),
               H[c, c'] = delta_cc' (1 + sigma/(2(1 - sigma - x_c'))) (c in U),  H[t, c'] = w_t w_c'/(2(1 - sigma - x_c'))
  Y_QQ = 2 (I + N'/2)^{-1},  N'[c,c] = sum_{u in U, u != c} x_u/(2 - sigma - 2x_c - 2x_u),
               N'[c,c'] = w_c w_c'/(2 - sigma - 2x_c - 2x_c')
Inverse-free reduction (congruence rho = (I + N1/2) rho', then elimination of q with Y_QQ^{-1} = (I + N'/2)/2):
  X11 = (I + N1/2) D_R (I + N1/2) - H (I + N'/2) H^T / 2
  X12 = (I + N1/2) Y_RP + H (I + N'/2) Y_PQ^T / 2
  X22 = Pi - Y_PQ (I + N'/2) Y_PQ^T / 2
The block is >= 0 with kernel exactly p (x) |2@t> iff X = [[X11, X12], [X12^T, X22]] is >= 0 with kernel exactly
span(0, w).  Checked: every closed-form block equals the Lemma B block; X equals the congruence-transformed Schur
complement of the Lemma B block; X (0, w) = 0.  Usage: python check_t1_closed.py"""
import random
from fractions import Fraction as Fr
from j1lib import J1, block_states, inv

random.seed(5)


def mat(n, m, f):
    return [[f(i, j) for j in range(m)] for i in range(n)]


def mul(A, B):
    return [[sum(A[i][l] * B[l][j] for l in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def T(A):
    return [list(r) for r in zip(*A)]


def add(A, B, c=1):
    return [[a + c * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def closed(w, t):
    k = len(w)
    x = [v * v for v in w]   # the caller ensures sum w^2 = 1
    U = [c for c in range(k) if c != t]
    sig = x[t]
    N = [1 if r == t else 2 for r in range(k)]
    g = [xc / (2 * (1 - xc)) for xc in x]
    gam = sum(g)
    tau = (1 + g[t]) / (1 + gam / 2)
    gamp = sum(x[u] / (2 - sig - 2 * x[u]) for u in U)
    taup = 1 / (1 + gamp / 2)
    W = w
    DR = mat(k, k, lambda i, j: Fr(N[i]) if i == j else Fr(0))
    taus = [tau - 2 if c == t else tau for c in range(k)]
    YRP = mat(k, k, lambda c, r: taus[c] / (2 * (1 - x[c])) * ((1 if r == c else 0) - W[r] * W[c]))
    Pi = mat(k, k, lambda r, rp: (N[r] if r == rp else 0) - N[r] * N[rp] * W[r] * W[rp] / (2 - sig))
    YPQ = mat(k, len(U), lambda r, j: taup / (2 - sig - 2 * x[U[j]]) * (N[r] * W[r] * W[U[j]] - (2 - sig) * (1 if r == U[j] else 0)))
    N1 = mat(k, k, lambda c, cp: sum(x[u] / (2 * (1 - x[c] - x[u])) for u in range(k) if u != c) if c == cp
             else W[c] * W[cp] / (2 * (1 - x[c] - x[cp])))
    H = mat(k, len(U), lambda c, j: ((1 + sig / (2 * (1 - sig - x[U[j]]))) if c == U[j] else 0) if c != t
            else W[t] * W[U[j]] / (2 * (1 - sig - x[U[j]])))
    I_k = mat(k, k, lambda i, j: Fr(1) if i == j else Fr(0))
    I_u = mat(len(U), len(U), lambda i, j: Fr(1) if i == j else Fr(0))
    A1 = add(I_k, N1, Fr(1, 2))
    Np = mat(len(U), len(U), lambda i, j: sum(x[u] / (2 - sig - 2 * x[U[i]] - 2 * x[u]) for u in U if u != U[i]) if i == j
             else W[U[i]] * W[U[j]] / (2 - sig - 2 * x[U[i]] - 2 * x[U[j]]))
    A2 = add(I_u, Np, Fr(1, 2))
    YRQ = [[-v for v in row] for row in mul(inv(A1, Fr), H)]
    YQQ = [[2 * v for v in row] for row in inv(A2, Fr)]
    return dict(DR=DR, YRP=YRP, Pi=Pi, YPQ=YPQ, YRQ=YRQ, YQQ=YQQ, A1=A1, A2=A2, H=H, U=U)


def unit_weights(k):
    """rational w with sum of squares 1: use w = 2v/(1 + |v|^2)-type parametrisation of the sphere."""
    v = [Fr(random.randint(1, 9), random.randint(1, 9)) for _ in range(k - 1)]
    n2 = sum(a * a for a in v)
    w = [2 * a / (1 + n2) for a in v] + [(1 - n2) / (1 + n2)]
    w = [abs(a) for a in w]
    if any(a == 0 for a in w):
        return unit_weights(k)
    return w


ok = True
cnt = 0
for k in range(3, 9):
    for trial in range(3):
        w = unit_weights(k)
        assert sum(a * a for a in w) == 1
        t = random.randrange(k)
        J = J1(w)
        st = block_states(k, (t,))
        M = J.matrix(st)
        Ridx = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 2]
        Pidx = [i for i, (a, y) in enumerate(st) if a is not None and a[0] == 1]
        Qidx = [i for i, (a, y) in enumerate(st) if a is None]
        # order: R_c by c, P_r by r, Q_c by c in U (block_states lists them in this order)
        assert [st[i][1].index(1) for i in Ridx] == list(range(k))
        Cf = closed(w, t)
        sub = lambda I, Jx: [[M[i][j] for j in Jx] for i in I]
        good = (sub(Ridx, Ridx) == Cf['DR'] and sub(Ridx, Pidx) == Cf['YRP'] and sub(Pidx, Pidx) == Cf['Pi']
                and sub(Pidx, Qidx) == Cf['YPQ'] and sub(Ridx, Qidx) == Cf['YRQ'] and sub(Qidx, Qidx) == Cf['YQQ'])
        # inverse-free reduction
        A1, A2, H = Cf['A1'], Cf['A2'], Cf['H']
        X11 = add(mul(mul(A1, Cf['DR']), A1), mul(mul(H, A2), T(H)), Fr(-1, 2))
        X12 = add(mul(A1, Cf['YRP']), mul(mul(H, A2), T(Cf['YPQ'])), Fr(1, 2))
        X22 = add(Cf['Pi'], mul(mul(Cf['YPQ'], A2), T(Cf['YPQ'])), Fr(-1, 2))
        # direct: congruence then Schur complement on q of the Lemma B block
        RP = Ridx + Pidx
        Bq = sub(RP, Qidx)
        Cq = inv(sub(Qidx, Qidx), Fr)
        Sch = add(sub(RP, RP), mul(mul(Bq, Cq), T(Bq)), -1)
        Dm = [[(A1[i][j] if (i < k and j < k) else (Fr(1) if i == j else Fr(0))) for j in range(2 * k)] for i in range(2 * k)]
        Xd = mul(mul(T(Dm), Sch), Dm)
        Xc = [X11[i] + X12[i] for i in range(k)] + [T(X12)[i] + X22[i] for i in range(k)]
        good2 = Xd == Xc
        kern = all(sum(Xc[i][k + r] * w[r] for r in range(k)) == 0 for i in range(2 * k))
        cnt += 1
        if not (good and good2 and kern):
            print(f'FAIL k={k} t={t}: closed form {good}, reduction {good2}, kernel {kern}')
        ok &= good and good2 and kern
    print(f'k = {k}: closed forms, inverse-free reduction and kernel exact: {ok}', flush=True)
print(f'{cnt} blocks checked')
print('PASS' if ok else 'FAIL')
