"""Code subspaces of the CSS code CSS(RM_3(2,2), RM_3(2,2)) on the 9 points of AG(2,3) (k = 3 logical qutrits).
For an encoded rank-two operator C = Enc(C_L):  q_9(t, C) = 3^{-9} sum_{logical Weyl L} W_{L+S}(1-3t) |Tr(C_L L)|^2,
W_X(lam) = sum_{E in X} lam^{9 - wt E}.  Minimise over rank-two C_L on (C^3)^{(x)3} (complex factors, L-BFGS).
usage: python rm_code.py t runs seed"""
import sys, itertools
import numpy as np
import torch
torch.set_default_dtype(torch.float64)
t, runs, seed = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
lam = 1 - 3 * t
n = 9
pts = [(x, y) for x in range(3) for y in range(3)]


def rm_gen(r):
    rows = []
    for a in range(3):
        for b in range(3):
            if a + b <= r:
                rows.append([(x ** a) * (y ** b) % 3 for (x, y) in pts])
    return np.array(rows, dtype=np.int64)


def rank3(M):
    M = M.copy() % 3; r = 0
    for c in range(M.shape[1]):
        piv = next((i for i in range(r, M.shape[0]) if M[i, c] % 3), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), -1, 3)) % 3
        for i in range(M.shape[0]):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % 3
        r += 1
    return r, M[:r]


G2 = rank3(rm_gen(2))[1]; G1 = rank3(rm_gen(1))[1]
print("dim RM(2) =", G2.shape[0], " dim RM(1) =", G1.shape[0], " RM(1) in RM(2):",
      rank3(np.vstack([G2, G1]))[0] == G2.shape[0], " RM(1) orthogonal to RM(2):", not ((G1 @ G2.T) % 3).any())
# logical X-bar_i, Z-bar_i: complete RM(1) to RM(2) by 3 vectors, then pair them by the dot product
comp = []
basis = G1.copy()
for v in G2:
    if rank3(np.vstack([basis, v]))[0] > basis.shape[0]:
        basis = np.vstack([basis, v]); comp.append(v)
comp = np.array(comp)
Gm = (comp @ comp.T) % 3                    # pairing matrix of the complement (mod RM(1))
Ginv = np.round(np.linalg.inv(Gm) * np.linalg.det(Gm)).astype(np.int64) * pow(int(round(np.linalg.det(Gm))) % 3, -1, 3) % 3
Zc = (Ginv.T @ comp) % 3                     # dual basis: Xc_i . Zc_j = delta_ij mod 3
assert ((comp @ Zc.T) % 3 == np.eye(3, dtype=np.int64)).all()


def span(gens, dimlen):
    P = np.zeros((1, dimlen), dtype=np.int64)
    for g in gens:
        P = np.concatenate([(P + k * g) % 3 for k in range(3)])
    return P


# stabilizer S = {X^a Z^b : a, b in RM(1)} ; vectors (x | z) in F_3^{18}
SX = span(G1, n); S = np.array([np.concatenate([a, b]) for a in SX for b in SX])
def wt(P):
    return ((P[:, :n] != 0) | (P[:, n:] != 0)).sum(1)
weightsN = None
w = {}
for a in itertools.product(range(3), repeat=3):
    for b in itertools.product(range(3), repeat=3):
        vx = (np.array(a) @ comp) % 3; vz = (np.array(b) @ Zc) % 3
        coset = (S + np.concatenate([vx, vz])) % 3
        w[(a, b)] = np.sum(lam ** (n - wt(coset)))
allN = np.concatenate([(S + np.concatenate([(np.array(a) @ comp) % 3, (np.array(b) @ Zc) % 3])) % 3
                       for a in itertools.product(range(3), repeat=3) for b in itertools.product(range(3), repeat=3)])
print("normalizer weight distribution:", np.bincount(wt(allN), minlength=n + 1).tolist())
print("stabilizer weight distribution:", np.bincount(wt(S), minlength=n + 1).tolist())
# logical Weyl operators on 3 qutrits
om = np.exp(2j * np.pi / 3)
X1 = np.roll(np.eye(3), 1, axis=0); Z1 = np.diag([1, om, om ** 2])
def weyl(a, b):
    M = np.array([[1.0 + 0j]])
    for i in range(3):
        M = np.kron(M, np.linalg.matrix_power(X1, a[i]) @ np.linalg.matrix_power(Z1, b[i]))
    return M
keys = list(w.keys())
Ls = torch.tensor(np.array([weyl(a, b) for (a, b) in keys]))        # 729 x 27 x 27
ws = torch.tensor(np.array([w[k] for k in keys]))


def qval(p):
    A = torch.complex(p[0], p[1]); B = torch.complex(p[2], p[3]); C = A @ B.conj().T
    tr = torch.einsum('kij,ji->k', Ls, C)
    return (ws * (tr.real ** 2 + tr.imag ** 2)).sum() / 3 ** 9 / (C.real ** 2 + C.imag ** 2).sum()


rng = np.random.default_rng(seed); best = np.inf
for r in range(runs):
    p = [torch.tensor(rng.standard_normal((27, 2)), requires_grad=True) for _ in range(4)]
    opt = torch.optim.LBFGS(p, lr=1, max_iter=2000, tolerance_grad=1e-13, tolerance_change=1e-16, history_size=50, line_search_fn='strong_wolfe')
    def cl():
        opt.zero_grad(); v = qval(p); v.backward(); return v
    for _ in range(3):
        opt.step(cl)
    best = min(best, qval(p).item())
print(f"t={t}: min over rank-two encoded C of q_9/||C||^2 = {best:+.6e}   (min logical weight {min(ws).item():+.4f}, identity weight {w[((0,0,0),(0,0,0))]:+.4f})")
