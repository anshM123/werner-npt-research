"""One free copy (copy 0) and one excitation on m sector copies: V = C^3 (x) (A + C e), A = span{0,1}^{(x)m}, e = 2...2.
 (1) Proposition M (exact identity with parameters rho, tau);
 (2) Lemma R': ran X in C^3 (x) c with rank X <= 2  =>  Sigma_full >= Gamma_X (adversarial, with a generic control);
 (3) Theorem E: ran C in C^3 (x) span{c, e} (c in A, possibly entangled), ANY co-range => q_n(C) >= 0 (adversarial, control);
 (4) the open inequality (*): sqrt(q1(Y)) (sqrt(Gamma_X) - sqrt(Sigma_full)) <= K sqrt(qt(P) qt(Q)); prints the largest
     ratio found (the sharp constant appears to be 4; the reduction needs K = 2^m - 2).
usage: python check_free_copy_excitation.py m starts"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
m, ST = int(sys.argv[1]), int(sys.argv[2]); torch.manual_seed(2026); n = m + 1; dA = 2 ** m; dK = dA + 1; D = 3 ** n
pa = [torch.eye(2, dtype=cplx), torch.tensor([[0, 1], [1, 0]], dtype=cplx), torch.tensor([[0, -1j], [1j, 0]], dtype=cplx), torch.diag(torch.tensor([1, -1], dtype=cplx))]
def kr(ms):
    out = torch.ones(1, 1, dtype=cplx)
    for M in ms: out = torch.kron(out, M)
    return out
SIG = {mu: torch.kron(torch.eye(3, dtype=cplx), kr([pa[s] for s in mu])) for mu in itertools.product(range(4), repeat=m)}
def q_dims(C, dims):
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0; letters = 'abcdefghijklmnop'
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
def q1(M): return torch.sum(torch.abs(M) ** 2) - 0.5 * torch.abs(torch.trace(M)) ** 2
cols = []
for a in range(3):
    for w in itertools.product(range(2), repeat=m):
        v = torch.zeros(D, dtype=cplx); v[int(np.ravel_multi_index((a,) + w, [3] * n))] = 1; cols.append(v)
    v = torch.zeros(D, dtype=cplx); v[int(np.ravel_multi_index((a,) + (2,) * m, [3] * n))] = 1; cols.append(v)
E = torch.stack(cols, 1)
def blocks(Cs):
    T = Cs.reshape(3, dK, 3, dK)
    return T[:, :dA, :, :dA].reshape(3 * dA, 3 * dA), T[:, dA, :, dA], T[:, :dA, :, dA], T[:, dA, :, :dA]
def xh(X, mu): return torch.einsum('iaja->ij', (SIG[mu] @ X).reshape(3, dA, 3, dA))
def sigma_full(X): return sum(q1(xh(X, mu)) for mu in itertools.product([1, 2, 3], repeat=m))
def gamma_X(X):
    T = X.reshape(3, dA, 3, dA)
    return 2 * sum(q1(T[:, a, :, b]) for a in range(dA) for b in range(dA)) - q1(torch.einsum('iaja->ij', T))
def qtP(P): return sum(q1(P[:, a, :]) for a in range(dA))
def qtQ(Q): return sum(q1(Q[:, :, a]) for a in range(dA))
# (1) Proposition M
err = 0.0
for _ in range(4):
    Cs = torch.randn(3 * dK, 3 * dK, dtype=cplx); rho, tau = float(np.exp(np.random.randn())), float(np.exp(np.random.randn()))
    X, Y, P, Q = blocks(Cs)
    D1 = torch.diag(torch.tensor([1.0] * dA + [(-1) ** (m + 1) * np.sqrt(rho / tau)], dtype=cplx)); D2 = torch.diag(torch.tensor([1.0] * dA + [np.sqrt(rho * tau)], dtype=cplx))
    Z = torch.kron(torch.eye(3, dtype=cplx), D1) @ Cs @ torch.kron(torch.eye(3, dtype=cplx), D2)
    lhs = 2 ** m * q_dims(E @ Cs @ E.conj().T, [3] * n) - (2 / rho) * q_dims(Z, [3, dK])
    rhs = sigma_full(X) - gamma_X(X) / rho + (1 - rho) * q1(Y) + (2 ** m - 2 * tau) * qtP(P) + (2 ** m - 2 / tau) * qtQ(Q)
    err = max(err, (abs(lhs - rhs) / torch.sum(torch.abs(Cs) ** 2)).item())
def adv(fn, npar, starts):
    best = 1e9
    for _ in range(starts):
        p = torch.randn(npar, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=600, tolerance_grad=1e-13, tolerance_change=1e-16, history_size=30, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = fn(p); v.backward(); return v
        try: opt.step(closure)
        except Exception: pass
        with torch.no_grad():
            v = fn(p).item()
            if np.isfinite(v): best = min(best, v)
    return best
def cvec(p, o, k): return torch.complex(p[o:o + k], p[o + k:o + 2 * k]), o + 2 * k
# (2) Lemma R'
def rp(p):
    c, o = cvec(p, 0, dA); Zl, o = cvec(p, o, 6); Zr, o = cvec(p, o, 6 * dA)
    X = torch.kron(torch.eye(3, dtype=cplx), c.reshape(dA, 1)) @ (Zl.reshape(3, 2) @ Zr.reshape(2, 3 * dA))
    return (sigma_full(X) - gamma_X(X)) / torch.sum(torch.abs(X) ** 2)
def rg(p):
    U, o = cvec(p, 0, 6 * dA); V, o = cvec(p, o, 6 * dA); X = U.reshape(3 * dA, 2) @ V.reshape(3 * dA, 2).conj().T
    return (sigma_full(X) - gamma_X(X)) / torch.sum(torch.abs(X) ** 2)
r_lemma = adv(rp, 2 * (dA + 6 + 6 * dA), ST); r_ctrl = adv(rg, 24 * dA, max(2, ST // 2))
# (3) Theorem E: range in C^3 (x) span{c, e}, arbitrary co-range
def te(p, rank=2):
    c, o = cvec(p, 0, dA); c = c / torch.linalg.norm(c)
    W = torch.zeros(3 * dK, 6, dtype=cplx)
    for a in range(3):
        W[a * dK:a * dK + dA, a] = c; W[a * dK + dA, 3 + a] = 1
    Ru, o = cvec(p, o, 6 * rank); Cv, o = cvec(p, o, D * rank)
    C = E @ W @ Ru.reshape(6, rank) @ Cv.reshape(D, rank).conj().T
    return q_dims(C, [3] * n) / torch.sum(torch.abs(C) ** 2)
r_E = adv(te, 2 * (dA + 12 + 2 * D), ST); r_E3 = adv(lambda p: te(p, 3), 2 * (dA + 18 + 3 * D), 1)
# (4) the open inequality (*): largest ratio over two-sided rank-2 members of the class
def star(p):
    U, o = cvec(p, 0, 6 * dK); V, o = cvec(p, o, 6 * dK); Cs = U.reshape(3 * dK, 2) @ V.reshape(3 * dK, 2).conj().T
    X, Y, P, Q = blocks(Cs)
    lhs = torch.sqrt(torch.clamp(q1(Y), min=1e-24)) * (torch.sqrt(torch.clamp(gamma_X(X), min=1e-24)) - torch.sqrt(torch.clamp(sigma_full(X), min=1e-24)))
    return -lhs / torch.sqrt(torch.clamp(qtP(P) * qtQ(Q), min=1e-24))
r_star = -adv(star, 24 * dK, ST)
print(f"m={m} (n={n}): Prop M error {err:.1e}; Lemma R' adversarial min gap {r_lemma:+.2e} (generic control {r_ctrl:+.3f}); "
      f"Theorem E adversarial min q_n/|C|^2 {r_E:+.2e} (rank-3 control {r_E3:+.3f}); largest ratio in (*) found {r_star:.4f} "
      f"(needed <= 2^m - 2 = {2 ** m - 2})")
