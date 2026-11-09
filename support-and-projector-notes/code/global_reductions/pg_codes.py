"""Two-dimensional subspaces of gauge-fixed projective Hamming codes (Section 7, after the RM_3(2,2) code).

PG(m-1,3) has N = (3^m - 1)/2 points; G (m x N) has the points as columns, generating the simplex code (every nonzero
word has weight 3^{m-1}).  Z-stabilizers: the simplex code.  X-stabilizers: X' = Hamming ∩ z^perp for a vector z not in
the simplex code, so X'^perp = <simplex, z>.  This is an [[N, 1, 3]]_3 code (k = 1) with logical Z = z and logical
X = x (G x = 0, z.x = 1).  By MacWilliams in the X part, the logical weights h_l = |S|^{-1} W_{l+S}(-1/2) are
    h(alpha, beta) = 3^{-m} sum_{b in beta z + simplex} 2^{|b|-N} sum_{d in <simplex, z>, d = 0 on supp b} w^{-alpha x.d} (-1)^{|d|},
and an encoded rank-two operator has q_N(C)/||C||_F^2 = sum_l h_l |Tr(C_L W_l^dag)|^2 / (3 ||C_L||_F^2).
The script checks the formula against a brute-force sum over all 3^12 stabilizers at N = 13, then reports, for random z,
the smallest q_N/||C||_F^2 over rank-two logical operators (local optimisation, floating point).
usage: python pg_codes.py m trials seed        (m = 3, 4, 5 give N = 13, 40, 121)"""
import sys, itertools
import numpy as np, torch
torch.set_default_dtype(torch.float64)
m, ntrials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(seed)
pts = [v for v in itertools.product(range(3), repeat=m) if any(v) and v[next(i for i in range(m) if v[i])] == 1]
G = np.array(pts, dtype=np.int64).T % 3
N = G.shape[1]
lam = -0.5; om = np.exp(2j * np.pi / 3)

def span(rows):
    rows = np.array(rows, dtype=np.int64) % 3
    co = np.array(list(itertools.product(range(3), repeat=rows.shape[0])), dtype=np.int64)
    return (co @ rows) % 3

def rref_solve(A, rhs):
    """solve A x = rhs over F_3 (one solution) or return None"""
    M = np.hstack([A % 3, (rhs % 3)[:, None]]); r = 0; piv = []
    for c in range(A.shape[1]):
        p = next((i for i in range(r, M.shape[0]) if M[i, c] % 3), None)
        if p is None: continue
        M[[r, p]] = M[[p, r]]; M[r] = (M[r] * (1 if M[r, c] == 1 else 2)) % 3
        for i in range(M.shape[0]):
            if i != r and M[i, c]: M[i] = (M[i] - M[i, c] * M[r]) % 3
        piv.append(c); r += 1
    if any(M[i, -1] for i in range(r, M.shape[0])): return None
    x = np.zeros(A.shape[1], dtype=np.int64)
    for i, c in enumerate(piv): x[c] = M[i, -1]
    return x

def nullspace(A):
    A = A.copy() % 3; r = 0; piv = []
    for c in range(A.shape[1]):
        p = next((i for i in range(r, A.shape[0]) if A[i, c] % 3), None)
        if p is None: continue
        A[[r, p]] = A[[p, r]]; A[r] = (A[r] * (1 if A[r, c] == 1 else 2)) % 3
        for i in range(A.shape[0]):
            if i != r and A[i, c]: A[i] = (A[i] - A[i, c] * A[r]) % 3
        piv.append(c); r += 1
    basis = []
    for f in [c for c in range(A.shape[1]) if c not in piv]:
        v = np.zeros(A.shape[1], dtype=np.int64); v[f] = 1
        for i, c in enumerate(piv): v[c] = (-A[i, f]) % 3
        basis.append(v)
    return np.array(basis)

def logical_x(z):
    return rref_solve(np.vstack([G, z]), np.concatenate([np.zeros(m, dtype=np.int64), [1]]))

def h_formula(z, x):
    Zp = span(G); Xperp = span(np.vstack([G, z]))
    sgn = (-1.0) ** (Xperp != 0).sum(1)
    h = np.zeros(9)
    for al in range(3):
        ph = om ** (-(al * (Xperp @ x)) % 3)
        for be in range(3):
            tot = 0
            for b in (Zp + be * z) % 3:
                sb = b != 0
                ok = ~(Xperp[:, sb] != 0).any(1)
                tot += 2.0 ** (sb.sum() - N) * (ph[ok] * sgn[ok]).sum()
            h[3 * al + be] = (tot / len(Zp)).real
    return h

def h_direct(z, x):
    gens = [np.concatenate([v, 0 * v]) for v in nullspace(np.vstack([G, z]))] + [np.concatenate([0 * g, g]) for g in G]
    S = span(gens)
    LX = np.concatenate([x, 0 * x]); LZ = np.concatenate([0 * z, z]); h = np.zeros(9)
    for al in range(3):
        for be in range(3):
            E = (S + al * LX + be * LZ) % 3
            h[3 * al + be] = (lam ** ((E[:, :N] == 0) & (E[:, N:] == 0)).sum(1)).sum() / len(S)
    return h

w = np.exp(2j * np.pi / 3)
Xm = np.roll(np.eye(3), 1, axis=0); Zm = np.diag([1, w, w * w])
Wt = torch.tensor(np.array([np.linalg.matrix_power(Xm, a) @ np.linalg.matrix_power(Zm, b)
                            for a in range(3) for b in range(3)]).conj())

def logical_min(h, restarts=16):
    ht = torch.tensor(h); best = None
    for _ in range(restarts):
        p = torch.tensor(rng.normal(size=24), requires_grad=True)
        def f(p):
            a = torch.complex(p[:6].reshape(2, 3), p[6:12].reshape(2, 3))
            b = torch.complex(p[12:18].reshape(2, 3), p[18:].reshape(2, 3))
            C = a.T @ b.conj(); c = torch.einsum('zij,ij->z', Wt, C)
            return (ht * c.abs() ** 2).sum() / (3 * (C.abs() ** 2).sum())
        opt = torch.optim.LBFGS([p], max_iter=200, line_search_fn='strong_wolfe')
        def cl():
            opt.zero_grad(); v = f(p); v.backward(); return v
        try: opt.step(cl)
        except Exception: pass
        v = f(p).item(); best = v if best is None else min(best, v)
    return best

if m == 3:
    z = rng.integers(0, 3, N); x = logical_x(z)
    while x is None:
        z = rng.integers(0, 3, N); x = logical_x(z)
    print(f"N = 13: formula versus brute force over 3^12 stabilizers: max difference "
          f"{np.abs(h_formula(z, x) - h_direct(z, x)).max():.1e}")
best = None
for trial in range(ntrials):
    while True:
        z = rng.integers(0, 3, N); x = logical_x(z)
        if x is not None: break
    h = h_formula(z, x)
    v = logical_min(h)
    best = v if best is None else min(best, v)
    print(f"N = {N}, trial {trial}: smallest weight {h.min():.3e}, largest {h.max():.3e}, "
          f"min q_N/||C||^2 over encoded rank-two C = {v:.4e}", flush=True)
print(f"N = {N}: smallest value {best:.4e} (positive; 2^(-N/3) = {2.0 ** (-N / 3):.3e})")
