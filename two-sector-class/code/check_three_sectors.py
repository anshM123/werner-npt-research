"""Three GHZ sectors: range/co-range in H_F (x) span{|1^m>, |2^m>, |3^m>} (random local bases at the m sector copies).
Claim: if m - g is even, q_n(C) = 2^{g-m} q_{F+{K_1..K_g}}(C_K) + (1 - 2^{g-m}) sum_{i!=j} q_F(X_ij),
with K_k = span{|1^{G_k}>, |2^{G_k}>, |3^{G_k}>} (a qutrit), C_K the same block matrix on H_F (x) span{|i>^{(x)g}}.
Then the adversarial minimum of q_n/|C|^2 over rank-two C in the class.
usage: python check_three_sectors.py f m g trials starts seed"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
f, m, g, TR, ST, SEED = [int(x) for x in sys.argv[1:7]]; n = f + m; DF = 3 ** f
torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
def q_dims(C, dims):
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0; letters = 'abcdefghijklmnopqrstuvwxyz'
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
Us = [torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]) for _ in range(m)]
def ghz_embed(bases, k):  # columns |alpha> (x) (x)_s bases[s][:, i], ordered (alpha, i)
    cols = []
    for a in range(DF):
        ea = torch.zeros(DF, dtype=cplx); ea[a] = 1
        for i in range(3):
            v = ea
            for s in range(k): v = torch.kron(v, bases[s][:, i])
            cols.append(v)
    return torch.stack(cols, 1)
E = ghz_embed(Us, m); EK = ghz_embed([torch.eye(3, dtype=cplx)] * g, g)
idmax = 0.0
for _ in range(TR):
    Cs = torch.randn(3 * DF, 3 * DF, dtype=cplx)
    lhs = q_dims(E @ Cs @ E.conj().T, [3] * n)
    blocks = Cs.reshape(DF, 3, DF, 3)
    off = sum(q_dims(blocks[:, i, :, j], [3] * f) if f > 0 else torch.abs(blocks[0, i, 0, j]) ** 2 for i in range(3) for j in range(3) if i != j)
    rhs = 2.0 ** (g - m) * q_dims(EK @ Cs @ EK.conj().T, [3] * (f + g)) + (1 - 2.0 ** (g - m)) * off
    idmax = max(idmax, (abs(lhs - rhs) / torch.sum(torch.abs(Cs) ** 2)).item())
best = 1e9; dV = 3 * DF
for _ in range(ST):
    p = torch.randn(8 * dV, requires_grad=True)
    def member(p):
        F_ = torch.complex(p[:2 * dV], p[2 * dV:4 * dV]).reshape(dV, 2); G_ = torch.complex(p[4 * dV:6 * dV], p[6 * dV:]).reshape(dV, 2)
        C = E @ F_ @ G_.conj().T @ E.conj().T; return C / torch.linalg.norm(C)
    opt = torch.optim.LBFGS([p], lr=1, max_iter=800, tolerance_grad=1e-13, tolerance_change=1e-17, history_size=40, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = q_dims(member(p), [3] * n); v.backward(); return v
    opt.step(closure)
    with torch.no_grad(): best = min(best, q_dims(member(p), [3] * n).item())
adv = f"; adversarial min q_n/|C|^2 over rank-2 members = {best:+.3e} ({ST} starts)" if ST > 0 else ""
print(f"f={f} m={m} g={g} (n={n}, m-g {'even' if (m-g)%2==0 else 'odd'}): identity max rel error {idmax:.1e}{adv}")
