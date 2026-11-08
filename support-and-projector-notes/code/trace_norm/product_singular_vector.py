"""Section 9, Proposition 26 (product singular vector): C = s1 |a><a| + s2 u v^*, a = (x) a_s product, u, v unit and orthogonal to a.
usage: python product_singular_vector.py n t runs
(i) check the identity G_x = sum_{K != 0} (-1)^|K| tau^(|K|-1) <v^(K), u^(K)>  (excitation sectors w.r.t. a);
(ii) adversarially minimise the TN gap over this class (n = 2, 3)."""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); t = float(sys.argv[2]); runs = int(sys.argv[3])
d = 3; D = d ** n
subsets = [K for r in range(n + 1) for K in itertools.combinations(range(n), r)]
def q(C):
    T = C.reshape([d] * 2 * n); tot = 0.
    for K in subsets:
        X, rem = T, list(range(n))
        for s in sorted(K, reverse=True):
            i = rem.index(s); X = torch.diagonal(X, dim1=i, dim2=len(rem) + i).sum(-1); rem.pop(i)
        tot = tot + (-t) ** len(K) * (X.real ** 2 + X.imag ** 2).sum()
    return tot
def build(p):
    k = 0
    def cvec(m):
        nonlocal k
        z = torch.complex(p[k], p[k + 1]); k += 2; return z / torch.linalg.norm(z)
    a_s = [cvec(d) for _ in range(n)]
    a = a_s[0]
    for x in a_s[1:]: a = torch.kron(a, x)
    u = torch.complex(p[k], p[k + 1]); v = torch.complex(p[k + 2], p[k + 3]); k += 4
    u = u - a * torch.vdot(a, u); v = v - a * torch.vdot(a, v)
    u = u / torch.linalg.norm(u); v = v / torch.linalg.norm(v)
    s1 = torch.exp(p[k][0]); s2 = torch.exp(p[k][1])
    return a_s, a, u, v, s1, s2
def gap(p):
    a_s, a, u, v, s1, s2 = build(p)
    C = s1 * torch.outer(a, a.conj()) + s2 * torch.outer(u, v.conj())
    n2 = (C.real ** 2 + C.imag ** 2).sum()
    rhs = (1 - t) ** (n - 1) * (n2 - t * (s1 + s2) ** 2)
    return (q(C) - rhs) / n2, (a_s, a, u, v, s1, s2)
rng = np.random.default_rng(5)
def init():
    ps = []
    for _ in range(n): ps += [torch.tensor(rng.standard_normal(d), requires_grad=True) for _ in range(2)]
    ps += [torch.tensor(rng.standard_normal(D), requires_grad=True) for _ in range(4)]
    ps += [torch.tensor(rng.standard_normal(2), requires_grad=True)]
    return ps
# (i) identity check on a random point
p = init()
with torch.no_grad():
    a_s, a, u, v, s1, s2 = build(p)
    tau = t / (1 - t)
    # G = tau^{-1}((I + 2 tau P_A)^{(x)n} - I) on H (x) H, compute G_x = <u1 v2|G|u2 v1> with u1 = v1 = a, u2 = u, v2 = v
    I9 = torch.eye(d * d, dtype=torch.complex128)
    F = torch.zeros(d * d, d * d, dtype=torch.complex128)
    for x in range(d):
        for y in range(d): F[x * d + y, y * d + x] = 1
    PA = (I9 - F) / 2
    Wn = torch.eye(1, dtype=torch.complex128)
    for s in range(n): Wn = torch.kron(Wn, I9 + 2 * tau * PA)
    # reorder from (A1 B1 A2 B2 ...) to (A1..An, B1..Bn)
    perm = [2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]
    Wn = Wn.reshape([d] * 4 * n).permute(perm + [2 * n + x for x in perm]).reshape(D * D, D * D)
    G = (Wn - torch.eye(D * D, dtype=torch.complex128)) / tau
    Gx = torch.vdot(torch.kron(a, v), G @ torch.kron(u, a))
    # excitation sectors w.r.t. a: basis at copy s = {a_s, two orthonormal vectors in a_s^perp}
    bases = []
    for s in range(n):
        M = torch.linalg.qr(torch.stack([a_s[s], torch.tensor(rng.standard_normal(d) + 0j), torch.tensor(rng.standard_normal(d) + 0j)], 1))[0]
        M = M * (torch.vdot(M[:, 0], a_s[s]) / abs(torch.vdot(M[:, 0], a_s[s])))  # first column = a_s up to phase -> fix phase
        bases.append(M)
    Bfull = bases[0]
    for M in bases[1:]: Bfull = torch.kron(Bfull, M)
    uc = Bfull.conj().T @ u; vc = Bfull.conj().T @ v      # coefficients in the adapted product basis
    tot = 0
    for idx in itertools.product(range(d), repeat=n):
        K = tuple(s for s in range(n) if idx[s] != 0)
        if not K: continue
        j = np.ravel_multi_index(idx, [d] * n)
        tot = tot + (-1) ** len(K) * tau ** (len(K) - 1) * vc[j].conj() * uc[j]
    print(f"(i) G_x = {Gx.real.item():+.10f} {Gx.imag.item():+.2e}i   excitation formula = {tot.real.item():+.10f} {tot.imag.item():+.2e}i")
best = 1e9
for r in range(runs):
    p = init()
    opt = torch.optim.LBFGS(p, lr=1, max_iter=2000, tolerance_grad=1e-12, tolerance_change=1e-15, history_size=50, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); g, _ = gap(p); g.backward(); return g
    for _ in range(3): opt.step(cl)
    best = min(best, gap(p)[0].item())
print(f"(ii) n={n} t={t}: smallest TN gap on the product-skeleton class = {best:.3e}")
