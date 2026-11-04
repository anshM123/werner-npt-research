"""Checks for Proposition 5 (the subsystem bound) of single-excitation-range/README.md.

Gauge sum_t A_t = I, sector c, h = n - c holes.  On the co-range sector V (x) H_c (index r*m + j, frame r, string j):
    Kt_c = I + Vac - T + 1/2 AA,   Vac = sum_t Abar_t (x) P_{0,t},  T = sum_t beta_t^* beta_t,  AA = alpha^* alpha,
    alpha = sum_t beta_t,  Ex = sum_t Abar_t (x) P_{1,t}  (= I - Vac stringwise),  L = 2T - AA.
(a) pair decomposition (exact):   L = 2/(h+1) sum_Q D_Q^* D_Q - (h-1)/(h+1) AA,
    Q over (c+1)-sets, (D_Q F)_{st} = beta_s P_{Q-t} F - beta_t P_{Q-s} F for s < t in Q (Theorem 3's operator on Q);
(b) Theorem 3 on a subsystem:     D_Q^* D_Q <= 2 Abar_Q (x) P_Q,   P_Q = projector onto strings with excited set in Q
    of size c (exactly one hole in Q);
(c) Proposition 5:                Kt_c >= (2h I - (3h-1) Ex)/(h+1) + (h-1)/(2(h+1)) AA,
    at random, structured and adversarial POVMs;
(d) the corollary: if lambda_max(A_S) <= 2h/(3h-1) for every c-set S then lambda_min Kt_c >= 0 (random POVMs that
    satisfy the condition).
"""
import numpy as np, itertools
from scipy.optimize import minimize

rng = np.random.default_rng(11)

def gauge(X):
    S = sum(x.conj().T @ x for x in X); w, V = np.linalg.eigh(S)
    return [x @ (V @ np.diag(w ** -0.5) @ V.conj().T) for x in X]

class Ops:
    def __init__(self, n, c):
        self.n, self.c = n, c
        lev = lambda k: [s for s in itertools.product(range(3), repeat=n) if sum(x != 0 for x in s) == k]
        self.S, self.R = lev(c), lev(c - 1)
        self.pS = {s: i for i, s in enumerate(self.S)}; self.pR = {s: i for i, s in enumerate(self.R)}
        self.m = len(self.S); self.supp = [frozenset(t for t in range(n) if s[t]) for s in self.S]

    def build(self, U):
        n, m = self.n, self.m; A = [u.conj().T @ u for u in U]; d = 2 * m
        beta = []
        for t in range(n):
            B = np.zeros((len(self.R), d), complex)
            for j, s in enumerate(self.S):
                if s[t] == 0: continue
                rho = list(s); rho[t] = 0; i = self.pR[tuple(rho)]
                for r in range(2): B[i, r * m + j] = np.conj(U[t][s[t] - 1, r])
            beta.append(B)
        def stringwise(fn):
            M = np.zeros((d, d), complex)
            for j in range(m):
                blk = fn(j)
                for r in range(2):
                    for rp in range(2): M[r * m + j, rp * m + j] = blk[r, rp]
            return M
        Vac = stringwise(lambda j: sum((A[t].T for t in range(n) if t not in self.supp[j]), np.zeros((2, 2))))
        Ex = stringwise(lambda j: sum((A[t].T for t in self.supp[j]), np.zeros((2, 2))))
        T = sum(b.conj().T @ b for b in beta); al = sum(beta); AA = al.conj().T @ al
        I = np.eye(d); K = I + Vac - T + 0.5 * AA
        return dict(A=A, beta=beta, Vac=Vac, Ex=Ex, T=T, AA=AA, I=I, K=(K + K.conj().T) / 2)

    def proj(self, X):   # projector onto strings whose excited set is X (on V (x) H_c)
        v = np.array([1.0 if self.supp[j] == X else 0.0 for j in range(self.m)]); return np.diag(np.concatenate([v, v]))

    def D(self, o, Q):
        rows = []
        for s, t in itertools.combinations(sorted(Q), 2):
            rows.append(o['beta'][s] @ self.proj(frozenset(Q) - {t}) - o['beta'][t] @ self.proj(frozenset(Q) - {s}))
        return np.vstack(rows)

    def AbarQ(self, o, Q):
        AQ = sum((o['A'][t].T for t in Q), np.zeros((2, 2))); m = self.m
        P = sum(self.proj(frozenset(Q) - {u}) for u in Q); M = np.zeros(P.shape, complex)
        for r in range(2):
            for rp in range(2): M[r * m:(r + 1) * m, rp * m:(rp + 1) * m] = AQ[r, rp] * P[:m, :m]
        return M

def prop5_gap(o, h):   # lambda_min of Kt_c - bound  (>= 0 required)
    B = (2 * h * o['I'] - (3 * h - 1) * o['Ex']) / (h + 1) + (h - 1) / (2 * (h + 1)) * o['AA']
    return np.linalg.eigvalsh(o['K'] - B)[0]

def random_U(n): return gauge([rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)])

ok = True
# (a), (b)
for n, c in [(3, 1), (4, 2), (4, 3), (5, 2), (5, 3)]:
    h = n - c; op = Ops(n, c); o = op.build(random_U(n))
    Qs = list(itertools.combinations(range(n), c + 1)); Ds = {Q: op.D(o, Q) for Q in Qs}
    lhs = 2 * o['T'] - o['AA']; rhs = 2 / (h + 1) * sum(D.conj().T @ D for D in Ds.values()) - (h - 1) / (h + 1) * o['AA']
    err = np.abs(lhs - rhs).max(); worst = max(np.linalg.eigvalsh(D.conj().T @ D - 2 * op.AbarQ(o, Q))[-1] for Q, D in Ds.items())
    print(f"n={n} c={c}: (a) pair decomposition error {err:.1e}; (b) max eigenvalue of D_Q^*D_Q - 2 Abar_Q P_Q {worst:+.1e}")
    ok &= err < 1e-12 and worst < 1e-12
# (c)
for n, c in [(3, 1), (3, 2), (4, 1), (4, 2), (4, 3), (5, 2), (5, 3), (5, 4), (6, 3), (6, 4)]:
    h = n - c; op = Ops(n, c); vals = []
    for k in range(6): vals.append(prop5_gap(op.build(random_U(n)), h))
    heavy = [np.array([[np.sqrt(0.9), 0], [0, 1.0]], complex)] + [np.array([[np.sqrt(0.1 / (n - 1)), 0], [0, 0]], complex)] * (n - 1)
    vals.append(prop5_gap(op.build(heavy), h)); vals.append(prop5_gap(op.build(gauge([np.eye(2, dtype=complex)] * n)), h))
    for k in range(2):
        X = [np.outer(rng.normal(size=2) + 1j * rng.normal(size=2), rng.normal(size=2)) + 0.03 * rng.normal(size=(2, 2)) for _ in range(n)]
        vals.append(prop5_gap(op.build(gauge(X)), h))
    if n <= 5:
        f = lambda x: prop5_gap(op.build(gauge([(x[4 * t:4 * t + 4] + 1j * x[4 * n + 4 * t:4 * n + 4 * t + 4]).reshape(2, 2) for t in range(n)])), h)
        for k in range(2):
            res = minimize(f, rng.normal(size=8 * n), method='L-BFGS-B', options={'maxiter': 60}); vals.append(res.fun)
    print(f"n={n} c={c} h={h}: (c) min over {len(vals)} POVMs of lambda_min(Kt_c - bound) = {min(vals):+.3e}")
    ok &= min(vals) > -1e-9
# (d)
for n, c in [(5, 3), (6, 3), (6, 4)]:
    h = n - c; op = Ops(n, c); thr = 2 * h / (3 * h - 1); found = 0; worst = np.inf; tries = 0
    while found < 6 and tries < 4000:
        tries += 1; eps = [0.3, 0.6, 1.0][tries % 3]       # near-uniform POVMs, so that the condition can hold
        U = gauge([np.eye(2) + eps * (rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))) for _ in range(n)]); A = [u.conj().T @ u for u in U]
        if max(np.linalg.eigvalsh(sum(A[t] for t in S))[-1] for S in itertools.combinations(range(n), c)) > thr: continue
        found += 1; worst = min(worst, np.linalg.eigvalsh(op.build(U)['K'])[0])
    print(f"n={n} c={c}: (d) {found} random POVMs with every c-set below {thr:.3f}; min lambda_min Kt_c = {worst:+.3e}")
    ok &= worst > -1e-12
assert ok
print("PASS")
