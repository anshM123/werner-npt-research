"""Products of GHZ_3 blocks: range/co-range in H_F (x) G_1 (x) ... (x) G_k, G_i = span{|1>^{B_i}, |2>^{B_i}, |3>^{B_i}} in random
local bases. Each block collapses to g_i qutrits (g_i = 1 for odd |B_i|, 2 for even). Check the two-stage identity (k = 2) and
the adversarial minimum of q_n/|C|^2 over rank-two members.  usage: python check_ghz_blocks.py f m1 m2 trials starts seed"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
f, m1, m2, TR, ST, SEED = [int(x) for x in sys.argv[1:7]]; torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
n = f + m1 + m2; DF = 3 ** f; g1 = 1 if m1 % 2 else 2; g2 = 1 if m2 % 2 else 2
def q_dims(C, dims):
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0; letters = 'abcdefghijklmnopqrstuvwxyz'
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
def haar(): return torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0])
def ghz_vecs(bases):  # 3 x (3^m) : row i = (x)_s bases[s][:, i]
    out = []
    for i in range(3):
        v = torch.ones(1, dtype=cplx)
        for Bs in bases: v = torch.kron(v, Bs[:, i])
        out.append(v)
    return torch.stack(out, 0)
G1 = ghz_vecs([haar() for _ in range(m1)]); G2 = ghz_vecs([haar() for _ in range(m2)])
K1 = ghz_vecs([torch.eye(3, dtype=cplx)] * g1); K2 = ghz_vecs([torch.eye(3, dtype=cplx)] * g2)
def embed(Ga, Gb):  # isometry C^{DF} (x) C^3 (x) C^3 -> H_F (x) span(Ga) (x) span(Gb)
    cols = []
    for a in range(DF):
        ea = torch.zeros(DF, dtype=cplx); ea[a] = 1
        for i in range(3):
            for j in range(3): cols.append(torch.kron(torch.kron(ea, Ga[i]), Gb[j]))
    return torch.stack(cols, 1)
E = embed(G1, G2); E_K1B2 = embed(K1, G2); E_K1K2 = embed(K1, K2)
def blocks1(Cs):  # sector blocks of block 1: operators on H_F (x) G_2 span, as 3 DF x 3 DF coefficient matrices
    T = Cs.reshape(DF, 3, 3, DF, 3, 3); return {(i, j): T[:, i, :, :, j, :].reshape(DF * 3, DF * 3) for i in range(3) for j in range(3)}
def blocks2(Cs):
    T = Cs.reshape(DF, 3, 3, DF, 3, 3); return {(i, j): T[:, :, i, :, :, j].reshape(DF * 3, DF * 3) for i in range(3) for j in range(3)}
EF_G2 = torch.stack([torch.kron(torch.eye(DF, dtype=cplx)[a], G2[j]) for a in range(DF) for j in range(3)], 1)
EF_K1 = torch.stack([torch.kron(torch.eye(DF, dtype=cplx)[a], K1[i]) for a in range(DF) for i in range(3)], 1)
idmax = 0.0
for _ in range(TR):
    Cs = torch.randn(9 * DF, 9 * DF, dtype=cplx)
    lhs = q_dims(E @ Cs @ E.conj().T, [3] * n)
    main = q_dims(E_K1K2 @ Cs @ E_K1K2.conj().T, [3] * (f + g1 + g2))
    b1 = blocks1(Cs); b2 = blocks2(Cs)
    off1 = sum(q_dims(EF_G2 @ b1[i, j] @ EF_G2.conj().T, [3] * (f + m2)) for i in range(3) for j in range(3) if i != j)
    off2 = sum(q_dims(EF_K1 @ b2[i, j] @ EF_K1.conj().T, [3] * (f + g1)) for i in range(3) for j in range(3) if i != j)
    c1 = 2.0 ** (g1 - m1); c2 = 2.0 ** (g2 - m2)
    rhs = c1 * c2 * main + c1 * (1 - c2) * off2 + (1 - c1) * off1
    idmax = max(idmax, (abs(lhs - rhs) / torch.sum(torch.abs(Cs) ** 2)).item())
best = 1e9; dV = 9 * DF
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
print(f"f={f} blocks=({m1},{m2}) -> g=({g1},{g2}), n={n}, collapsed copies {f + g1 + g2}: two-stage identity max rel error {idmax:.1e}; "
      f"adversarial min q_n/|C|^2 over rank-2 = {best:+.3e} ({ST} starts)")
