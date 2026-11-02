"""Theorem E with two free copies: F = copies 0,1 (qutrits), one excitation on m sector copies (A = span{0,1}^m, e = 2..2).
Checks: Proposition M with q_F (two-copy form on F) and q_{F+K}; Lemma R' with q_F; adversarial min of q_n over rank-2 C
with range in C^9 (x) span{c, e} (c in A random/entangled), any co-range.   usage: python check_free_copies_excitation.py m starts"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
m, ST = int(sys.argv[1]), int(sys.argv[2]); torch.manual_seed(77); n = m + 2; dA = 2 ** m; dK = dA + 1; D = 3 ** n; DF = 9
pa = [torch.eye(2, dtype=cplx), torch.tensor([[0, 1], [1, 0]], dtype=cplx), torch.tensor([[0, -1j], [1j, 0]], dtype=cplx), torch.diag(torch.tensor([1, -1], dtype=cplx))]
def kr(ms):
    out = torch.ones(1, 1, dtype=cplx)
    for M in ms: out = torch.kron(out, M)
    return out
SIG = {mu: torch.kron(torch.eye(DF, dtype=cplx), kr([pa[s] for s in mu])) for mu in itertools.product(range(4), repeat=m)}
def q_dims(C, dims):
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0; letters = 'abcdefghijklmnop'
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
def qF(M): return q_dims(M, [3, 3])
cols = []
for a in range(DF):
    a0, a1 = divmod(a, 3)
    for w in itertools.product(range(2), repeat=m):
        v = torch.zeros(D, dtype=cplx); v[int(np.ravel_multi_index((a0, a1) + w, [3] * n))] = 1; cols.append(v)
    v = torch.zeros(D, dtype=cplx); v[int(np.ravel_multi_index((a0, a1) + (2,) * m, [3] * n))] = 1; cols.append(v)
E = torch.stack(cols, 1)
def blocks(Cs):
    T = Cs.reshape(DF, dK, DF, dK)
    return T[:, :dA, :, :dA].reshape(DF * dA, DF * dA), T[:, dA, :, dA], T[:, :dA, :, dA], T[:, dA, :, :dA]
def xh(X, mu): return torch.einsum('iaja->ij', (SIG[mu] @ X).reshape(DF, dA, DF, dA))
def sigma_full(X): return sum(qF(xh(X, mu)) for mu in itertools.product([1, 2, 3], repeat=m))
def gamma_X(X):
    T = X.reshape(DF, dA, DF, dA)
    return 2 * sum(qF(T[:, a, :, b]) for a in range(dA) for b in range(dA)) - qF(torch.einsum('iaja->ij', T))
def qtP(P): return sum(qF(P[:, a, :]) for a in range(dA))
def qtQ(Q): return sum(qF(Q[:, :, a]) for a in range(dA))
err = 0.0
for _ in range(2):
    Cs = torch.randn(DF * dK, DF * dK, dtype=cplx); rho, tau = float(np.exp(np.random.randn())), float(np.exp(np.random.randn()))
    X, Y, P, Q = blocks(Cs)
    D1 = torch.diag(torch.tensor([1.0] * dA + [(-1) ** (m + 1) * np.sqrt(rho / tau)], dtype=cplx)); D2 = torch.diag(torch.tensor([1.0] * dA + [np.sqrt(rho * tau)], dtype=cplx))
    Z = torch.kron(torch.eye(DF, dtype=cplx), D1) @ Cs @ torch.kron(torch.eye(DF, dtype=cplx), D2)
    lhs = 2 ** m * q_dims(E @ Cs @ E.conj().T, [3] * n) - (2 / rho) * q_dims(Z, [3, 3, dK])
    rhs = sigma_full(X) - gamma_X(X) / rho + (1 - rho) * qF(Y) + (2 ** m - 2 * tau) * qtP(P) + (2 ** m - 2 / tau) * qtQ(Q)
    err = max(err, (abs(lhs - rhs) / torch.sum(torch.abs(Cs) ** 2)).item())
def adv(fn, npar, starts):
    best = 1e9
    for _ in range(starts):
        p = torch.randn(npar, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=400, tolerance_grad=1e-13, tolerance_change=1e-16, history_size=30, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = fn(p); v.backward(); return v
        try: opt.step(closure)
        except Exception: pass
        with torch.no_grad():
            v = fn(p).item()
            if np.isfinite(v): best = min(best, v)
    return best
def cvec(p, o, k): return torch.complex(p[o:o + k], p[o + k:o + 2 * k]), o + 2 * k
def rp(p):
    c, o = cvec(p, 0, dA); Zl, o = cvec(p, o, 2 * DF); Zr, o = cvec(p, o, 2 * DF * dA)
    X = torch.kron(torch.eye(DF, dtype=cplx), c.reshape(dA, 1)) @ (Zl.reshape(DF, 2) @ Zr.reshape(2, DF * dA))
    return (sigma_full(X) - gamma_X(X)) / torch.sum(torch.abs(X) ** 2)
def te(p, rank=2):
    c, o = cvec(p, 0, dA); c = c / torch.linalg.norm(c)
    W = torch.zeros(DF * dK, 2 * DF, dtype=cplx)
    for a in range(DF):
        W[a * dK:a * dK + dA, a] = c; W[a * dK + dA, DF + a] = 1
    Ru, o = cvec(p, o, 2 * DF * rank); Cv, o = cvec(p, o, D * rank)
    C = E @ W @ Ru.reshape(2 * DF, rank) @ Cv.reshape(D, rank).conj().T
    return q_dims(C, [3] * n) / torch.sum(torch.abs(C) ** 2)
r_lemma = adv(rp, 2 * (dA + 2 * DF + 2 * DF * dA), ST)
r_E = adv(te, 2 * (dA + 4 * DF + 2 * D), ST); r_E3 = adv(lambda p: te(p, 3), 2 * (dA + 6 * DF + 3 * D), 1)
print(f"two free copies, m={m} (n={n}): Prop M (q_F) error {err:.1e}; Lemma R' (q_F) adversarial min gap {r_lemma:+.2e}; "
      f"Theorem E2 adversarial min q_n/|C|^2 {r_E:+.2e} (rank-3 control {r_E3:+.3f})")
