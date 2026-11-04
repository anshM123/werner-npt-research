"""Structured balanced certificates for Proposition 7 (single-excitation-range/README.md).
Balanced certificate: matrices k_{S,u} (S a set of c copies, u a hole of S) with
  (P) sum_u k_{S,u} <= R_S = 2I + Abar_H - T_S          (every pattern S),
  (C) Pi_{S'} + (+)_{u in S'} k_{S'-u,u} >= 0               (every child S').
Here we restrict to three numbers per pair: k_{S,u} = x Abar_u / Tr A_u + w E_S + y I, with E_S = 2 Abar_S - T_S.
(a) n = 5, c = 3 and n = 6, c = 3 at structured real POVMs, and the heavy-copy family (restricted model) n = 5..10:
    solve for the numbers (maximise the slack of (C)), rebuild the matrices, re-check (P) and (C) with numpy, and
    repair solver tolerance by k -> k - eps I (eps = |min eig of (P)| / h + 1e-12), re-checking again;
(b) a closed-form rule (trace-weighted split of 2I + 2Abar_H - T_S) fails at some of these POVMs.
Optional: python check_structured_split.py 6 4   (n = 6, c = 4; slow, uses SCS)."""
import numpy as np, sys, itertools, warnings
import cvxpy as cp
from threelevel_certificate import GeneralModel, Model, gauge
warnings.filterwarnings('ignore')

def setup(M, c):
    o = M.ops(c); S_list = o['S']; m = o['m']; n = M.n
    supp = lambda s: frozenset(t for t in range(n) if s[t])
    pats = sorted(set(supp(s) for s in S_list), key=sorted)
    idx = {P: [r * m + j for r in range(2) for j, s in enumerate(S_list) if supp(s) == P] for P in pats}
    Abar = []
    for u in range(n):
        A = np.zeros((2 * m, 2 * m))
        for j in range(m):
            for r in range(2):
                for rp in range(2): A[r * m + j, rp * m + j] = np.real(M.A[u].conj()[r, rp])
        Abar.append(A)
    Et = []
    for t in range(n):
        exc = np.diag([1.0 if S_list[j][t] else 0.0 for r in range(2) for j in range(m)])
        Et.append(np.real(2 * exc @ Abar[t] @ exc - o['B'][t].conj().T @ o['B'][t]))
    R2 = np.real(2 * np.eye(2 * m) + o['Vac'] - o['T'])
    ch = []
    for Sp, Pi in o['Pi'].items():
        mem = [u for u in sorted(Sp) if (Sp - {u}) in idx]; ind = sum((idx[Sp - {u}] for u in mem), [])
        ch.append((Sp, mem, ind, np.real(Pi[np.ix_(ind, ind)])))
    return o, pats, idx, Abar, Et, R2, ch

def build(M, c, pre, V):
    o, pats, idx, Abar, Et, R2, ch = pre; n = M.n; blk = lambda X, P: X[np.ix_(idx[P], idx[P])]; K = {}
    for P in pats:
        ES = sum(blk(Et[t], P) for t in P)
        for u in [u for u in range(n) if u not in P]:
            x, w, y = V[P, u]; K[P, u] = x * blk(Abar[u], P) / max(np.trace(M.A[u]).real, 1e-12) + w * ES + y * np.eye(len(idx[P]))
    return K

def check(M, c, pre, K):
    o, pats, idx, Abar, Et, R2, ch = pre; n = M.n; blk = lambda X, P: X[np.ix_(idx[P], idx[P])]
    pmin = min(np.linalg.eigvalsh(blk(R2, P) - sum(K[P, u] for u in range(n) if u not in P))[0] for P in pats)
    cmin = np.inf
    for Sp, mem, ind, Pis in ch:
        Mx = Pis.copy(); off = 0
        for u in mem:
            d = len(idx[Sp - {u}]); Mx[off:off + d, off:off + d] += K[Sp - {u}, u]; off += d
        cmin = min(cmin, np.linalg.eigvalsh((Mx + Mx.T) / 2)[0])
    return pmin, cmin

def solve(M, c, pre):
    o, pats, idx, Abar, Et, R2, ch = pre; n = M.n; blk = lambda X, P: X[np.ix_(idx[P], idx[P])]
    lam = cp.Variable(); cons = []; K = {}; V = {}
    for P in pats:
        H = [u for u in range(n) if u not in P]; ES = sum(blk(Et[t], P) for t in P)
        for u in H:
            V[P, u] = (cp.Variable(), cp.Variable(), cp.Variable()); x, w, y = V[P, u]
            K[P, u] = x * blk(Abar[u], P) / max(np.trace(M.A[u]).real, 1e-12) + w * ES + y * np.eye(len(idx[P]))
        cons.append(blk(R2, P) - sum(K[P, u] for u in H) >> 0)
    for Sp, mem, ind, Pis in ch:
        rows = [[K[Sp - {u}, u] if i == j else np.zeros((len(idx[Sp - {u}]), len(idx[Sp - {mem[j]}]))) for j in range(len(mem))]
                for i, u in enumerate(mem)]
        cons.append(Pis + cp.bmat(rows) - lam * np.eye(len(ind)) >> 0)
    prob = cp.Problem(cp.Maximize(lam), cons)
    big = max(len(ind) for Sp, mem, ind, Pis in ch) > 100
    try:
        if big: raise RuntimeError
        prob.solve(solver='CLARABEL')
    except Exception:
        prob.solve(solver='SCS', eps_abs=1e-7, eps_rel=1e-7, max_iters=40000)
    return lam.value, {k: tuple(float(z.value) for z in v) for k, v in V.items()}

def run(name, M, c):
    pre = setup(M, c); lam, V = solve(M, c, pre); K = build(M, c, pre, V); p0, c0 = check(M, c, pre, K)
    h = M.n - c; eps = max(0.0, -p0) / h + 1e-12
    K = {k: v - eps * np.eye(v.shape[0]) for k, v in K.items()}; p1, c1 = check(M, c, pre, K)
    lm = 2 * np.linalg.eigvalsh(M.ops(c)['K'])[0]
    print(f"(a) n={M.n} c={c} {name:22s} solver slack {lam:+.3e}; re-check (P) {p0:+.1e} -> {p1:+.1e} after eps = {eps:.1e}; "
          f"(C) {c1:+.3e}; lambda_min(2Kt) {lm:+.3e} >= h (C) = {h * c1:+.3e}", flush=True)
    return p1 >= 0 and c1 > 0 and lm >= h * c1 - 1e-9

def cases(n, rng):
    def rz(*s): return rng.normal(size=s)
    out = [("random real", gauge([rz(2, 2) for _ in range(n)]))]
    X = [np.outer(rz(2), rz(2)) + 0.05 * rz(2, 2) for _ in range(n)]; out.append(("near rank-one frame", gauge(X)))
    X = [np.outer(rz(2), rz(2)) + 0.05 * rz(2, 2) for _ in range(n)]; X[0] *= 0.02; X[1] *= 0.02; out.append(("two nearly unused copies", gauge(X)))
    X = [np.eye(2) + 0.1 * rz(2, 2), np.eye(2) + 0.1 * rz(2, 2)] + [0.1 * rz(2, 2) for _ in range(n - 2)]; out.append(("two heavy copies", gauge(X)))
    X = [np.diag([1.0, 1.0]) + 0.03 * rz(2, 2)] + [0.15 * rz(2, 2) for _ in range(n - 1)]; out.append(("one heavy copy", gauge(X)))
    out.append(("uniform W-plane", gauge([np.eye(2)] * n)))
    return out

if __name__ == '__main__':
    ok = True
    if len(sys.argv) > 2:
        n, c = int(sys.argv[1]), int(sys.argv[2])
        for name, U in cases(n, np.random.default_rng(606)): ok &= run(name, GeneralModel(U), c)
    else:
        for name, U in cases(5, np.random.default_rng(7)): ok &= run(name, GeneralModel(U), 3)
        for name, U in cases(6, np.random.default_rng(606)): ok &= run(name, GeneralModel(U), 3)
        for N, d in [(4, 0.1), (6, 0.1), (6, 0.2), (7, 0.2), (8, 0.2), (9, 0.2)]: ok &= run(f"heavy family d={d}", Model(N, d), N - 1)
        # (b) closed-form rule: G_{S,u} = (TrA_u / TrA_H) (2I + 2Abar_H - T_S), i.e. k = G - Abar_u
        worst = np.inf; rng = np.random.default_rng(606)
        for name, U in cases(6, rng)[:3]:
            M = GeneralModel(U); pre = setup(M, 4); o, pats, idx, Abar, Et, R2, ch = pre; n = M.n; K = {}
            blk = lambda X, P: X[np.ix_(idx[P], idx[P])]; a = [np.trace(x).real for x in M.A]
            for P in pats:
                H = [u for u in range(n) if u not in P]; DS = blk(R2, P) + sum(blk(Abar[u], P) for u in H)
                for u in H: K[P, u] = a[u] / sum(a[x] for x in H) * DS - blk(Abar[u], P)
            worst = min(worst, check(M, 4, pre, K)[1])
        print(f"(b) trace-weighted closed-form split at n = 6, c = 4 (three POVMs): smallest child eigenvalue {worst:+.3e} < 0")
        ok &= worst < 0
    assert ok
    print("PASS")
