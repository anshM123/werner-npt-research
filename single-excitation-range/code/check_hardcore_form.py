"""Checks for the remarks after Proposition 5 (single-excitation-range/README.md).

Gauge sum_t A_t = I, co-range sector c, F = (G^1, G^2) in V (x) H_c.  a_r^dag = sum_t |w_{t,r}><0|_t creates the range
vector u_r = sum_t (U_t e_r)@t, a_r is its adjoint, alpha F = a_1 G^1 + a_2 G^2.
(i)   hard-core form:  <F|Kt_c|F> = ||F||^2 + sum_{r,r'} <a_{r'}^dag G^r, a_r^dag G^{r'}> - 1/2 ||alpha F||^2;
(ii)  twisted identity (no gauge needed):  L = (sum_t Tr A_t) ||F||^2 - slack - sum_t ||gamma_t Ft||^2 - ||Theta_+||^2
      + 1/2 ||sum_t gamma_t Ft||^2,  Ft = (G^2, -G^1), gamma_t = sum_r b_{t,r}^dag <e_r|,
      slack = sum_t [Tr A_t ||P_{1,t} F||^2 - ||beta_t F||^2],
      ||Theta_+||^2 = ||a_1^dag G^1||^2 + ||a_2^dag G^2||^2 + 1/2 ||a_1^dag G^2 + a_2^dag G^1||^2;
(iii) c = 1: every created pair has two terms, so 1/2 ||sum_t gamma_t Ft||^2 <= sum_t ||gamma_t Ft||^2 and L <= 2||F||^2;
(iv)  the vacuum terms are needed: one heavy copy A_0 = diag(1-d, 1) plus n-1 light copies (d/(n-1)) |e_1><e_1| gives
      L > 2||F||^2 for some F in every sector 2 <= c <= n-2 (n = 5, 6), while Kt_c >= 0.
"""
import numpy as np, itertools
rng = np.random.default_rng(5)

def gauge(X):
    S = sum(x.conj().T @ x for x in X); w, V = np.linalg.eigh(S)
    return [x @ (V @ np.diag(w ** -0.5) @ V.conj().T) for x in X]

def op_at(M, t, n):
    out = np.eye(1)
    for s in range(n): out = np.kron(out, M if s == t else np.eye(3))
    return out

def ops(U, n):
    b = {}; P0 = {}; P1 = {}
    for t in range(n):
        for r in range(2):
            w = np.zeros(3, complex); w[1:] = U[t][:, r]; b[t, r] = op_at(np.outer(np.eye(3)[0], w.conj()), t, n)
        P0[t] = op_at(np.diag([1.0, 0, 0]), t, n); P1[t] = op_at(np.diag([0, 1.0, 1.0]), t, n)
    return b, P0, P1

def level(n, c): return np.array([1.0 if sum(x != 0 for x in s) == c else 0.0 for s in itertools.product(range(3), repeat=n)])
def n2(x): return np.vdot(x, x).real

ok = True; worst_i = worst_ii = 0
for n in (3, 4):
    U = gauge([rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]); A = [u.conj().T @ u for u in U]
    b, P0, P1 = ops(U, n)
    ad = [sum(b[t, r].conj().T for t in range(n)) for r in range(2)]
    for c in range(1, n + 1):
        lp = level(n, c)
        for trial in range(3):
            G = [(rng.normal(size=3 ** n) + 1j * rng.normal(size=3 ** n)) * lp for _ in range(2)]; nF = n2(G[0]) + n2(G[1])
            beta = [b[t, 0] @ G[0] + b[t, 1] @ G[1] for t in range(n)]; T = sum(n2(x) for x in beta); al = sum(beta)
            vac = sum(np.vdot(G[r], P0[t] @ G[rp]) * np.conj(A[t][r, rp]) for t in range(n) for r in range(2) for rp in range(2)).real
            K = nF + vac - T + 0.5 * n2(al)                                   # <F|Kt_c|F>
            ThS = sum(np.vdot(ad[rp] @ G[r], ad[r] @ G[rp]) for r in range(2) for rp in range(2)).real
            worst_i = max(worst_i, abs(K - (nF + ThS - 0.5 * n2(al))) / nF)
            L = 2 * T - n2(al); Gt = [G[1], -G[0]]
            gt = [b[t, 0].conj().T @ Gt[0] + b[t, 1].conj().T @ Gt[1] for t in range(n)]
            slack = sum(np.trace(A[t]).real * (n2(P1[t] @ G[0]) + n2(P1[t] @ G[1])) - n2(beta[t]) for t in range(n))
            thp = n2(ad[0] @ G[0]) + n2(ad[1] @ G[1]) + 0.5 * n2(ad[0] @ G[1] + ad[1] @ G[0])
            rhs = sum(np.trace(a).real for a in A) * nF - slack - sum(n2(x) for x in gt) - thp + 0.5 * n2(sum(gt))
            worst_ii = max(worst_ii, abs(L - rhs) / nF)
            if c == 1: ok &= 0.5 * n2(sum(gt)) <= sum(n2(x) for x in gt) + 1e-9 and L <= 2 * nF + 1e-9
print(f"(i) hard-core form, n = 3, 4, all c: largest relative error {worst_i:.1e}")
print(f"(ii) twisted identity, n = 3, 4, all c: largest relative error {worst_ii:.1e}")
print(f"(iii) c = 1 pair bound and L <= 2||F||^2 at all samples: {ok}")
ok &= worst_i < 1e-12 and worst_ii < 1e-12

# (iv) heavy copy: sector operator without and with the vacuum terms (explicit sector formula)
def sector(U, n, c, with_vac):
    strs = [s for s in itertools.product(range(3), repeat=n) if sum(x != 0 for x in s) == c]; pos = {s: i for i, s in enumerate(strs)}; m = len(strs)
    A = [u.conj().T @ u for u in U]; K = np.zeros((2 * m, 2 * m), complex)
    for j, s in enumerate(strs):
        vac = [t for t in range(n) if s[t] == 0]
        K[j::m, j::m] += np.eye(2) + (sum((A[t].T for t in vac), np.zeros((2, 2))) if with_vac else 0)
        for t in range(n):
            if s[t] == 0: continue
            for zp in (1, 2):
                sp = list(s); sp[t] = zp; K[pos[tuple(sp)]::m, j::m] += -0.5 * np.outer(U[t][zp - 1, :], U[t][s[t] - 1, :].conj())
            for u in vac:
                for z in (1, 2):
                    sp = list(s); sp[t] = 0; sp[u] = z; K[pos[tuple(sp)]::m, j::m] += 0.5 * np.outer(U[u][z - 1, :], U[t][s[t] - 1, :].conj())
    return np.linalg.eigvalsh((K + K.conj().T) / 2)[0]
for n in (5, 6):
    for c in range(2, n - 1):
        best = None
        for d in np.linspace(0.02, 0.5, 13):
            U = [np.array([[np.sqrt(1 - d), 0], [0, 1.0]], complex)] + [np.array([[np.sqrt(d / (n - 1)), 0], [0, 0]], complex)] * (n - 1)
            nv, wv = sector(U, n, c, False), sector(U, n, c, True)
            if best is None or nv < best[0]: best = (nv, wv, d)
        print(f"(iv) n={n} c={c}: heavy copy d={best[2]:.2f}: without vacuum lambda_min {best[0]:+.3e} (< 0: L > 2||F||^2), "
              f"with vacuum {best[1]:+.3e}")
        ok &= best[0] < 0 and best[1] > 0
assert ok
print("PASS")
