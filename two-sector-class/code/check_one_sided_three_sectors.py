"""One-sided three-sector identity: ran C in V = H_F (x) span{1^m, 2^m, 3^m} (random local bases), co-range ARBITRARY.
Claim: q_n(C) = q_n(P_V C P_V) + sum_{i, w non-constant} 2^{-|A(i,w)|} q_F(<i^m| C |w>),  A(i,w) = {s : w_s = i}.
Then adversarial min of q_n/|C|^2 over rank-2 C with range in V (any co-range), and of the products G3 x G3 one-sided.
usage: python check_one_sided_three_sectors.py f m trials starts seed"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
f, m, TR, ST, SEED = [int(x) for x in sys.argv[1:6]]; n = f + m; DF = 3 ** f; D = 3 ** n
torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
def q_dims(C, dims):
    k = len(dims)
    if k == 0: return torch.abs(C.reshape(())) ** 2 if C.numel() == 1 else torch.sum(torch.abs(C) ** 2)
    T = C.reshape(dims + dims); tot = 0.0; letters = 'abcdefghijklmnopqrstuvwxyz'
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
Us = [torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]) for _ in range(m)]
def string(w):
    v = torch.ones(1, dtype=cplx)
    for s, c in enumerate(w): v = torch.kron(v, Us[s][:, c])
    return v
Ufull = torch.ones(1, 1, dtype=cplx)
for U in Us: Ufull = torch.kron(Ufull, U)
Ufull = torch.kron(torch.eye(DF, dtype=cplx), Ufull)            # columns: |alpha> (x) |w> in local bases, w in [3]^m
words = list(itertools.product(range(3), repeat=m)); widx = {w: k for k, w in enumerate(words)}
const = [widx[(i,) * m] for i in range(3)]
cols = [a * 3 ** m + const[i] for a in range(DF) for i in range(3)]
EV = Ufull[:, cols]; PV = EV @ EV.conj().T
idmax = 0.0
for _ in range(TR):
    Cs = torch.randn(3 * DF, D, dtype=cplx)                       # range coefficients (alpha, i) x arbitrary co-range
    C = EV @ Cs                                                  # D x D operator with range in V
    lhs = q_dims(C, [3] * n); main = q_dims(PV @ C @ PV, [3] * n)
    Cl = C @ Ufull                                               # co-range expressed in the local basis (columns |alpha'> (x) |w>)
    extra = 0.0
    for i in range(3):
        for w in words:
            if len(set(w)) == 1: continue
            Aw = sum(1 for c in w if c == i)
            rows = [a * 3 ** m + const[i] for a in range(DF)]; cl = [a * 3 ** m + widx[w] for a in range(DF)]
            Mf = (Ufull[:, rows].conj().T @ C @ Ufull[:, cl])     # <i^m| C |w> as an operator on H_F
            extra = extra + 2.0 ** (-Aw) * q_dims(Mf, [3] * f)
    idmax = max(idmax, (abs(lhs - main - extra) / torch.linalg.norm(C) ** 2).item())
best = 1e9; dR = 3 * DF
for _ in range(ST):
    p = torch.randn(4 * dR + 4 * D, requires_grad=True)
    def member(p):
        F_ = torch.complex(p[:2 * dR], p[2 * dR:4 * dR]).reshape(dR, 2); G_ = torch.complex(p[4 * dR:4 * dR + 2 * D], p[4 * dR + 2 * D:]).reshape(D, 2)
        C = EV @ F_ @ G_.conj().T; return C / torch.linalg.norm(C)
    opt = torch.optim.LBFGS([p], lr=1, max_iter=800, tolerance_grad=1e-13, tolerance_change=1e-17, history_size=40, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = q_dims(member(p), [3] * n); v.backward(); return v
    opt.step(closure)
    with torch.no_grad(): best = min(best, q_dims(member(p), [3] * n).item())
print(f"f={f} m={m} (n={n}): one-sided identity max rel error {idmax:.1e}; adversarial min q_n/|C|^2 (range in V, any co-range) = {best:+.3e} ({ST} starts)")
