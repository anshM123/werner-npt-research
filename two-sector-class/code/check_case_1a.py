"""One free copy + one excitation, case (1a): range and co-range each contain a nonzero vector of C^3 (x) A.
Checks Lemma T on random two-piece X (sqrt Gamma - sqrt Sigma <= sqrt Gamma(E) + sqrt Sigma(E)), the bounds Sigma <= 2^m |Z|^2 and
Gamma <= 2|Z|^2 (rank <= 2), the (*) ratio in case (1a) against the proved constant 2 sqrt2 + 2^{1+m/2}, and the adversarial
minimum of q_n over case-(1a) members.   usage: python check_case_1a.py m starts"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
m, ST = int(sys.argv[1]), int(sys.argv[2]); torch.manual_seed(13); n = m + 1; dA = 2 ** m; dK = dA + 1; D = 3 ** n
pa = [torch.eye(2, dtype=cplx), torch.tensor([[0, 1], [1, 0]], dtype=cplx), torch.tensor([[0, -1j], [1j, 0]], dtype=cplx), torch.diag(torch.tensor([1, -1], dtype=cplx))]
def kr(ms):
    out = torch.ones(1, 1, dtype=cplx)
    for M in ms: out = torch.kron(out, M)
    return out
SIG = [torch.kron(torch.eye(3, dtype=cplx), kr([pa[s] for s in mu])) for mu in itertools.product(range(1, 4), repeat=m)]
def q1(M): return torch.sum(torch.abs(M) ** 2) - 0.5 * torch.abs(torch.trace(M)) ** 2
def sig(X): return sum(q1(torch.einsum('iaja->ij', (S @ X).reshape(3, dA, 3, dA))) for S in SIG)
def gam(X):
    T = X.reshape(3, dA, 3, dA)
    return 2 * sum(q1(T[:, a, :, b]) for a in range(dA) for b in range(dA)) - q1(torch.einsum('iaja->ij', T))
def q_dims(C, dims):
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0; letters = 'abcdefghijklmnop'
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
cols = []
for a in range(3):
    for w in itertools.product(range(2), repeat=m):
        v = torch.zeros(D, dtype=cplx); v[int(np.ravel_multi_index((a,) + w, [3] * n))] = 1; cols.append(v)
    v = torch.zeros(D, dtype=cplx); v[int(np.ravel_multi_index((a,) + (2,) * m, [3] * n))] = 1; cols.append(v)
E = torch.stack(cols, 1)
worstT, worstS, worstG = -1e9, -1e9, -1e9
for _ in range(20):
    x1, x2, y1, y2 = [torch.randn(3 * dA, dtype=cplx) for _ in range(4)]
    Z1 = torch.outer(x1, y1.conj()); X0 = torch.outer(x2, y2.conj()); X = X0 + Z1
    lhs = torch.sqrt(torch.clamp(gam(X), min=0)) - torch.sqrt(torch.clamp(sig(X), min=0))
    rhs = torch.sqrt(torch.clamp(gam(Z1), min=0)) + torch.sqrt(torch.clamp(sig(Z1), min=0))
    worstT = max(worstT, (lhs - rhs).item())
    worstS = max(worstS, (sig(X) / torch.sum(torch.abs(X) ** 2)).item() / 2 ** m); worstG = max(worstG, (gam(X) / torch.sum(torch.abs(X) ** 2)).item() / 2)
def blocks(Cs):
    T = Cs.reshape(3, dK, 3, dK)
    return T[:, :dA, :, :dA].reshape(3 * dA, 3 * dA), T[:, dA, :, dA], T[:, :dA, :, dA], T[:, dA, :, :dA]
def member(p):   # case (1a): U_B = u e1^T, V_B = v e1^T
    o = 0
    def take(k):
        nonlocal o
        z = torch.complex(p[o:o + k], p[o + k:o + 2 * k]); o += 2 * k; return z
    UA = take(6 * dA).reshape(3 * dA, 2); VA = take(6 * dA).reshape(3 * dA, 2); u = take(3); v = take(3)
    UB = torch.stack([u, torch.zeros(3, dtype=cplx)], 1); VB = torch.stack([v, torch.zeros(3, dtype=cplx)], 1)
    U = torch.cat([UA.reshape(3, dA, 2), UB.reshape(3, 1, 2)], 1).reshape(3 * dK, 2); V = torch.cat([VA.reshape(3, dA, 2), VB.reshape(3, 1, 2)], 1).reshape(3 * dK, 2)
    return U @ V.conj().T
npar = 2 * (12 * dA + 6)
def star_ratio(p):
    X, Y, P, Q = blocks(member(p))
    qp = sum(q1(P[:, a, :].reshape(3, 3)) for a in range(dA)); qq = sum(q1(Q[:, :, a].reshape(3, 3)) for a in range(dA))
    lhs = torch.sqrt(torch.clamp(q1(Y), min=1e-24)) * (torch.sqrt(torch.clamp(gam(X), min=1e-24)) - torch.sqrt(torch.clamp(sig(X), min=1e-24)))
    return -lhs / torch.sqrt(torch.clamp(qp * qq, min=1e-24))
def qn(p):
    C = member(p); return q_dims(E @ C @ E.conj().T, [3] * n) / torch.sum(torch.abs(C) ** 2)
def adv(fn, starts):
    best = 1e9
    for _ in range(starts):
        p = torch.randn(npar, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=0.5, max_iter=400, tolerance_grad=1e-12, tolerance_change=1e-15, history_size=30, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = fn(p); v.backward(); return v
        try: opt.step(closure)
        except Exception: pass
        with torch.no_grad():
            v = fn(p).item()
            if np.isfinite(v): best = min(best, v)
    return best
rmax = -adv(star_ratio, ST); qmin = adv(qn, max(2, ST // 2))
print(f"m={m}: Lemma T max violation {worstT:+.1e}; max Sigma/(2^m|X|^2) {worstS:.3f}; max Gamma/(2|X|^2) {worstG:.3f} (random two-piece X); "
      f"case (1a): largest (*) ratio {rmax:.4f} vs proved constant {2 * np.sqrt(2) + 2 ** (1 + m / 2):.2f} and needed {2 ** m - 2}; adversarial min q_n/|C|^2 {qmin:+.2e}")
