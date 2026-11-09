"""Section 7 (projective Hamming codes), without gauge fixing. Full [[13,7,3]]_3 code: CSS(Hamming ⊇ simplex) on PG(2,3).  All 3^14 logical weights h(alpha, beta) at lambda = -1/2:
  h = 3^-3 sum_{b in beta.z + simplex} 2^{|b|-13} sum_{d in Hamming, d = 0 on supp b} w^{-alpha.(x.d)} (-1)^{|d|}
computed with a subset-sum (zeta) transform over supports and a ternary Fourier transform in alpha.
Then (if some weight is negative) minimise q/||C||^2 = sum_L h_L |Tr(W_L^dag C)|^2 / (3^7 ||C||^2) over rank-two C on the
2187-dimensional logical space.  A negative value is a counterexample to EP_13.
usage: python hamming13_full.py restarts seed"""
import sys, itertools, time
import numpy as np, torch
torch.set_default_dtype(torch.float64)
restarts, seed = int(sys.argv[1]), int(sys.argv[2])
rng = np.random.default_rng(seed); torch.manual_seed(seed)
m = 3
pts = [v for v in itertools.product(range(3), repeat=m) if any(v) and v[next(i for i in range(m) if v[i])] == 1]
G = np.array(pts, dtype=np.int64).T % 3; N = G.shape[1]; k = N - 2 * m           # 13, 7
def span(rows):
    rows = np.array(rows, dtype=np.int64) % 3
    co = np.array(list(itertools.product(range(3), repeat=rows.shape[0])), dtype=np.int64)
    return (co @ rows) % 3
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
def inv_mod3(M):
    n = M.shape[0]; A = np.hstack([M % 3, np.eye(n, dtype=np.int64)])
    for c in range(n):
        p = next(i for i in range(c, n) if A[i, c] % 3)
        A[[c, p]] = A[[p, c]]; A[c] = (A[c] * (1 if A[c, c] == 1 else 2)) % 3
        for i in range(n):
            if i != c and A[i, c]: A[i] = (A[i] - A[i, c] * A[c]) % 3
    return A[:, n:]
simplex = span(G)
Hb = nullspace(G)                                        # 10 x 13 basis of the Hamming code
# logical basis: Hamming vectors independent modulo the simplex code
def rank_mod3(A):
    A = np.array(A, dtype=np.int64) % 3; r = 0
    for c in range(A.shape[1]):
        p = next((i for i in range(r, A.shape[0]) if A[i, c] % 3), None)
        if p is None: continue
        A[[r, p]] = A[[p, r]]; A[r] = (A[r] * (1 if A[r, c] == 1 else 2)) % 3
        for i in range(A.shape[0]):
            if i != r and A[i, c]: A[i] = (A[i] - A[i, c] * A[r]) % 3
        r += 1
    return r
U = []
for v in Hb:
    if rank_mod3(np.vstack(list(G) + U + [v])) == m + len(U) + 1:
        U.append(v)
    if len(U) == k: break
assert len(U) == k, "could not find a complement of the simplex code"
U = np.array(U); Mg = (U @ U.T) % 3
try:
    Minv = inv_mod3(Mg)
except StopIteration:
    raise SystemExit("logical Gram matrix singular mod 3; pick another basis")
X = U; Z = (Minv.T @ U) % 3
assert (((X @ Z.T) % 3) == np.eye(k, dtype=np.int64)).all() and not ((X @ G.T) % 3).any()
Hall = span(Hb)                                           # 3^10 Hamming codewords
supp = ((Hall != 0) * (1 << np.arange(N))).sum(1)
sgn = np.where((Hall != 0).sum(1) % 2 == 0, 1, -1).astype(np.int32)
iota = (((Hall @ X.T) % 3) * (3 ** np.arange(k))).sum(1)
t0 = time.time()
F = np.zeros((1 << N, 3 ** k), dtype=np.int32)
np.add.at(F, (supp, iota), sgn)
Fz = F.reshape([2] * N + [3 ** k])                        # bit i of the support index = axis (N-1-i)
for ax in range(N):
    Fz = np.cumsum(Fz, axis=ax, dtype=np.int32)          # subset sums: sum over s ⊆ U
Gz = Fz.reshape(1 << N, 3 ** k)
print(f"zeta transform done ({time.time() - t0:.0f}s)", flush=True)
om = np.exp(2j * np.pi / 3)
D3 = np.array([[om ** (-(a * b) % 3) for b in range(3)] for a in range(3)])        # alpha x iota digit
def ft_alpha(vec):
    T = vec.reshape([3] * k)                              # axis order: digit k-1 ... 0 (C order of reshape)
    for ax in range(k):
        T = np.moveaxis(np.tensordot(D3, T, axes=([1], [ax])), 0, ax)
    return T.reshape(-1)
h = np.zeros((3 ** k, 3 ** k))                            # h[alpha_index, beta_index]
full = (1 << N) - 1
for bi, beta in enumerate(itertools.product(range(3), repeat=k)):
    beta = np.array(beta[::-1])                           # beta index digits little-endian
    zvec = (beta @ Z) % 3
    acc = np.zeros(3 ** k)
    for s in simplex:
        b = (zvec + s) % 3; mask = int(((b != 0) * (1 << np.arange(N))).sum())
        acc += 2.0 ** (bin(mask).count('1') - N) * Gz[full ^ mask]
    h[:, bi] = (ft_alpha(acc) / len(simplex)).real
print(f"weights done ({time.time() - t0:.0f}s): min {h.min():.4e}, max {h.max():.4e}, "
      f"negative entries {int((h < -1e-15).sum())} of {h.size}", flush=True)
# sanity: identity weight directly
Sfull = [(sx, sz) for sx in simplex for sz in simplex]
direct = np.mean([(-0.5) ** int(((sx == 0) & (sz == 0)).sum()) for sx, sz in Sfull])
print(f"check: identity weight formula {h[0, 0]:.6e} vs direct {direct:.6e}")
alpha_t = torch.tensor(np.array(list(itertools.product(range(3), repeat=k)))[:, ::-1].copy())

if (h < -1e-15).sum() == 0:
    print("all logical weights are nonnegative: the encoded form is positive semidefinite, no counterexample here")
    raise SystemExit
D = 3 ** k
digits = np.array(list(itertools.product(range(3), repeat=k)))[:, ::-1]          # index -> little-endian digits
plus = np.zeros((D, D), dtype=np.int64)
for a in range(D):
    plus[a] = ((((digits + digits[a]) % 3) * (3 ** np.arange(k))).sum(1))
plus_t = torch.tensor(plus)
ht = torch.tensor(h)
D3t = torch.tensor(D3)
def coeffs(C):
    V = torch.gather(C, 0, plus_t.T)                      # V[a, j] = C[j + a, j]
    T = V.reshape([D] + [3] * k)                          # (alpha, digits of j)
    for ax in range(1, k + 1):
        T = torch.movedim(torch.tensordot(D3t, T, dims=([1], [ax])), 0, ax)
    return T.reshape(D, D)                                # [alpha, beta] = Tr(W^dag C)
# checks: Parseval, the identity, and rank-one operators (q > 0 for rank one, by the trace-norm bound)
Cr = torch.complex(torch.tensor(rng.normal(size=(D, D))), torch.tensor(rng.normal(size=(D, D))))
cr = coeffs(Cr)
print(f"Parseval: sum|c|^2 / (D ||C||^2) = {((cr.abs() ** 2).sum() / (D * (Cr.abs() ** 2).sum())).item():.12f}")
Ci = torch.eye(D, dtype=torch.complex128)
print(f"identity: value {((ht * coeffs(Ci).abs() ** 2).sum() / (D * D)).item():+.6e} (should equal h[0,0] = {h[0, 0]:+.6e})")
vals1 = []
for _ in range(5):
    a = torch.complex(torch.tensor(rng.normal(size=D)), torch.tensor(rng.normal(size=D)))
    C1 = torch.outer(a, a.conj())
    vals1.append(((ht * coeffs(C1).abs() ** 2).sum() / (D * (C1.abs() ** 2).sum())).item())
print(f"random rank-one values (must be > 0): {['%+.3e' % v for v in vals1]}")
best = None
for r in range(restarts):
    p = torch.tensor(rng.normal(size=8 * D) / np.sqrt(D), requires_grad=True)
    def f(p):
        a = torch.complex(p[:2 * D].reshape(2, D), p[2 * D:4 * D].reshape(2, D))
        b = torch.complex(p[4 * D:6 * D].reshape(2, D), p[6 * D:].reshape(2, D))
        C = a.T @ b.conj()
        c = coeffs(C)
        return (ht * c.abs() ** 2).sum() / (D * (C.abs() ** 2).sum())
    opt = torch.optim.LBFGS([p], max_iter=150, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); v = f(p); v.backward(); return v
    opt.step(cl)
    v = f(p).item(); best = v if best is None else min(best, v)
    print(f"restart {r}: q_13/||C||^2 = {v:+.6e} ({time.time() - t0:.0f}s)", flush=True)
print(f"smallest q_13/||C||^2 over encoded rank-two operators: {best:+.6e}  "
      f"({'NEGATIVE -> COUNTEREXAMPLE CANDIDATE' if best < -1e-10 else 'no counterexample'})")
