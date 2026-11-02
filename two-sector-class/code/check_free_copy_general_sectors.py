"""Theorem E' (one free copy, general two-sector structure with arbitrary a, b): range in C^3 (x) span{c, c'}, c in A, c' in B,
any co-range => q_n >= 0.  Checks the generalised Proposition M
  2^m q_n(C) - (2/rho) q_{0,K}(Z) = (Sigma_X - Gamma_X/rho) + (Sigma_Y - rho Gamma_Y) + (2^m - 2 tau) qt(P) + (2^m - 2/tau) qt(Q)
and the adversarial minimum.   usage: python check_free_copy_general_sectors.py types starts   (types like 21,12,21)"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
TYPES = [tuple(int(c) for c in t) for t in sys.argv[1].split(',')]; ST = int(sys.argv[2]); torch.manual_seed(5); rng = np.random.default_rng(5)
m = len(TYPES); n = m + 1; D = 3 ** n
Us = [torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]) for _ in range(m)]
def kron_all(ms):
    out = torch.ones(1, 1, dtype=cplx)
    for M in ms: out = torch.kron(out, M)
    return out
SA = kron_all([Us[s][:, :TYPES[s][0]] for s in range(m)]); SB = kron_all([Us[s][:, TYPES[s][0]:sum(TYPES[s])] for s in range(m)])
dA, dB = SA.shape[1], SB.shape[1]; dK = dA + dB
E = torch.kron(torch.eye(3, dtype=cplx), torch.cat([SA, SB], 1))        # C^3 (x) (A + B) -> (C^3)^n
pa = [torch.eye(2, dtype=cplx), torch.tensor([[0, 1], [1, 0]], dtype=cplx), torch.tensor([[0, -1j], [1j, 0]], dtype=cplx), torch.diag(torch.tensor([1, -1], dtype=cplx))]
def paulis(side):   # Pauli strings acting on the planar factors of A (side 0) or B (side 1), identity on lines
    planar = [s for s in range(m) if TYPES[s][side] == 2]
    out = {}
    for mu in itertools.product(range(4), repeat=len(planar)):
        ms = []; it = iter(mu)
        for s in range(m):
            dim = TYPES[s][side]
            ms.append(pa[next(it)] if dim == 2 else torch.eye(1, dtype=cplx))
        out[mu] = torch.kron(torch.eye(3, dtype=cplx), kron_all(ms))
    return out
PA, PB = paulis(0), paulis(1)
def q_dims(C, dims):
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0; letters = 'abcdefghijklmnop'
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
def q1(M): return torch.sum(torch.abs(M) ** 2) - 0.5 * torch.abs(torch.trace(M)) ** 2
def hat(Xs, S, d): return torch.einsum('iaja->ij', (S @ Xs).reshape(3, d, 3, d))
def sig(Xs, PS, d): return sum(q1(hat(Xs, S, d)) for mu, S in PS.items() if all(x > 0 for x in mu))
def gam(Xs, d):
    T = Xs.reshape(3, d, 3, d)
    return 2 * sum(q1(T[:, a, :, b]) for a in range(d) for b in range(d)) - q1(torch.einsum('iaja->ij', T))
def blocks(Cs):
    T = Cs.reshape(3, dK, 3, dK)
    return (T[:, :dA, :, :dA].reshape(3 * dA, 3 * dA), T[:, dA:, :, dA:].reshape(3 * dB, 3 * dB), T[:, :dA, :, dA:], T[:, dA:, :, :dA])
def qtPQ(B4, d1, d2): return sum(q1(B4[:, a, :, b]) for a in range(d1) for b in range(d2))
err = 0.0
for _ in range(3):
    Cs = torch.randn(3 * dK, 3 * dK, dtype=cplx); rho, tau = float(np.exp(np.random.randn())), float(np.exp(np.random.randn()))
    X, Y, P, Q = blocks(Cs); eps = (-1) ** (m + 1)
    D1 = torch.diag(torch.tensor([1.0] * dA + [eps * np.sqrt(rho / tau)] * dB, dtype=cplx)); D2 = torch.diag(torch.tensor([1.0] * dA + [np.sqrt(rho * tau)] * dB, dtype=cplx))
    Z = torch.kron(torch.eye(3, dtype=cplx), D1) @ Cs @ torch.kron(torch.eye(3, dtype=cplx), D2)
    lhs = 2 ** m * q_dims(E @ Cs @ E.conj().T, [3] * n) - (2 / rho) * q_dims(Z, [3, dK])
    rhs = (sig(X, PA, dA) - gam(X, dA) / rho) + (sig(Y, PB, dB) - rho * gam(Y, dB)) + (2 ** m - 2 * tau) * qtPQ(P, dA, dB) + (2 ** m - 2 / tau) * qtPQ(Q, dB, dA)
    err = max(err, (abs(lhs - rhs) / torch.sum(torch.abs(Cs) ** 2)).item())
def cvec(p, o, k): return torch.complex(p[o:o + k], p[o + k:o + 2 * k]), o + 2 * k
def te(p, rank=2):
    c, o = cvec(p, 0, dA); c2, o = cvec(p, o, dB); c = c / torch.linalg.norm(c); c2 = c2 / torch.linalg.norm(c2)
    W = torch.zeros(3 * dK, 6, dtype=cplx)
    for a in range(3):
        W[a * dK:a * dK + dA, a] = c; W[a * dK + dA:a * dK + dK, 3 + a] = c2
    Ru, o = cvec(p, o, 6 * rank); Cv, o = cvec(p, o, D * rank)
    C = E @ W @ Ru.reshape(6, rank) @ Cv.reshape(D, rank).conj().T
    return q_dims(C, [3] * n) / torch.sum(torch.abs(C) ** 2)
def adv(fn, npar, starts):
    best = 1e9
    for _ in range(starts):
        p = torch.randn(npar, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=500, tolerance_grad=1e-13, tolerance_change=1e-16, history_size=30, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = fn(p); v.backward(); return v
        try: opt.step(closure)
        except Exception: pass
        with torch.no_grad():
            v = fn(p).item()
            if np.isfinite(v): best = min(best, v)
    return best
r = adv(te, 2 * (dA + dB + 12 + 2 * D), ST); r3 = adv(lambda p: te(p, 3), 2 * (dA + dB + 18 + 3 * D), 1)
print(f"types {TYPES} (n={n}, a={sum(1 for t in TYPES if t[0]==2)}, b={sum(1 for t in TYPES if t[1]==2)}): generalised Prop M error {err:.1e}; "
      f"Theorem E' adversarial min q_n/|C|^2 {r:+.2e} ({ST} starts; rank-3 control {r3:+.3f})")
