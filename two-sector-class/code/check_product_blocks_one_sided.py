"""Adversarial test of the block-collapse theorem on mixed products.  Blocks (comma separated):
  F      one free qutrit copy (weight 1)
  G<m>   GHZ_3 block span{1^m, 2^m, 3^m} in random local bases (weight 1 if m odd, 2 if even)
  T<ab>  two-sector block A + B with A = (x) A_s, B = (x) B_s, types like 21.12.11 (one group: <= 1 plane per side)
ONE-SIDED: range of C in the product of the blocks' subspaces, co-range arbitrary.  Prints min q_n/|C|^2 over rank two.
usage: python check_product_blocks_one_sided.py BLOCKS starts seed      e.g.  T21.12,G3   or  F,G3,G1"""
import sys, numpy as np, torch
torch.set_default_dtype(torch.float64); cplx = torch.complex128
spec, ST, SEED = sys.argv[1].split(','), int(sys.argv[2]), int(sys.argv[3]); torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
def q_dims(C, dims):
    k = len(dims); T = C.reshape(dims + dims); tot = 0.0; letters = 'abcdefghijklmnopqrstuvwxyz'
    for mask in range(2 ** k):
        J = [i for i in range(k) if mask >> i & 1]; keep = [i for i in range(k) if i not in J]
        L = list(letters[:k]); R = [letters[k + i] if i in keep else letters[i] for i in range(k)]
        out = ''.join(L[i] for i in keep) + ''.join(R[i] for i in keep)
        tot = tot + (-0.5) ** len(J) * torch.sum(torch.abs(torch.einsum(''.join(L) + ''.join(R) + '->' + out, T)) ** 2)
    return tot
def haar(): return torch.tensor(np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0])
def kron_all(vs):
    v = torch.ones(1, dtype=cplx)
    for w in vs: v = torch.kron(v, w)
    return v
blocks, ncop, weight = [], 0, 0
for b in spec:
    if b == 'F':
        blocks.append(torch.eye(3, dtype=cplx)); ncop += 1; weight += 1
    elif b[0] == 'G':
        m = int(b[1:]); Us = [haar() for _ in range(m)]
        blocks.append(torch.stack([kron_all([U[:, i] for U in Us]) for i in range(3)], 1)); ncop += m; weight += 1 if m % 2 else 2
    elif b[0] == 'T':
        types = [(int(t[0]), int(t[1])) for t in b[1:].split('.')]; Us = [haar() for _ in types]
        def sector(side):
            vecs = [torch.ones(1, dtype=cplx)]
            for (da, db), U in zip(types, Us):
                cols = U[:, :da] if side == 0 else U[:, da:da + db]
                vecs = [torch.kron(v, cols[:, c]) for v in vecs for c in range(cols.shape[1])]
            return vecs
        blocks.append(torch.stack(sector(0) + sector(1), 1)); ncop += len(types)
        weight += max(sum(1 for t in types if t[0] == 2), sum(1 for t in types if t[1] == 2), 1)
E = blocks[0]
for Bk in blocks[1:]: E = torch.kron(E, Bk)
dV = E.shape[1]
best = 1e9
for _ in range(ST):
    DD = E.shape[0]; p = torch.randn(4 * dV + 4 * DD, requires_grad=True)
    def member(p):
        F_ = torch.complex(p[:2 * dV], p[2 * dV:4 * dV]).reshape(dV, 2); G_ = torch.complex(p[4 * dV:4 * dV + 2 * DD], p[4 * dV + 2 * DD:]).reshape(DD, 2)
        C = E @ F_ @ G_.conj().T; return C / torch.linalg.norm(C)
    opt = torch.optim.LBFGS([p], lr=1, max_iter=800, tolerance_grad=1e-13, tolerance_change=1e-17, history_size=40, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = q_dims(member(p), [3] * ncop); v.backward(); return v
    opt.step(closure)
    with torch.no_grad(): best = min(best, q_dims(member(p), [3] * ncop).item())
# control: the optimiser does find negatives for rank three on the same support
DD = E.shape[0]; p3 = torch.randn(6 * dV + 6 * DD, requires_grad=True)
def member3(p):
    F_ = torch.complex(p[:3 * dV], p[3 * dV:6 * dV]).reshape(dV, 3); G_ = torch.complex(p[6 * dV:6 * dV + 3 * DD], p[6 * dV + 3 * DD:]).reshape(DD, 3)
    C = E @ F_ @ G_.conj().T; return C / torch.linalg.norm(C)
opt = torch.optim.LBFGS([p3], lr=1, max_iter=800, tolerance_grad=1e-13, tolerance_change=1e-17, history_size=40, line_search_fn='strong_wolfe')
def closure3():
    opt.zero_grad(); v = q_dims(member3(p3), [3] * ncop); v.backward(); return v
opt.step(closure3)
with torch.no_grad(): r3 = q_dims(member3(p3), [3] * ncop).item()
print(f"blocks {spec}: n={ncop}, support dim {dV}, weight sum {weight}: adversarial min q_n/|C|^2 over rank 2 = {best:+.3e} ({ST} starts); rank-3 control {r3:+.3f}")
