"""Checks for the one-excitation class (README.md).  K_s = span{|0>,|1>} at every copy, e = |2...2>, V = K^n + C e.
(1) the identity  q_n(C) = q_n(C_K) + |kappa|^2 + |kappa'|^2 + 2^-n |g|^2 + 2 (-1/2)^n Re(conj(g) Tr C_K)
    on random rank-two C with range and co-range in V;
(2) the qubit lemma  2^n q_n(X) >= (|Tr X| - 2^n s_2(X))_+^2  for rank-two X on (C^2)^n, adversarially;
(3) the minimum of q_n(C)/||C||^2 over the class (rank two), adversarially;
(4) the two-copy thick critical point code/critical_n2.npy lies in such a class (after a local unitary).
usage: python check_one_excitation.py [nmax]"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
cplx = torch.complex128
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
torch.manual_seed(0); np.random.seed(0)
w9 = torch.eye(3, dtype=cplx).reshape(9) / np.sqrt(3)
PAULI = [torch.tensor(m, dtype=cplx) for m in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])]
def make_tools(n):
    D = 3 ** n
    def Phi(C):
        X = C.reshape((3,) * (2 * n)).permute([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n)
        for s in range(n):
            tr = torch.tensordot(w9.conj(), X, dims=([0], [s])); X = X - 1.5 * torch.moveaxis(torch.tensordot(w9, tr, dims=0), 0, s)
        return X.reshape([3, 3] * n).permute([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]).reshape(D, D)
    def q(C): return torch.vdot(C.reshape(-1), Phi(C).reshape(-1)).real
    Kidx = [int(''.join(map(str, t)), 3) for t in itertools.product([0, 1], repeat=n)]; eidx = D - 1
    fullP = []
    for mus in itertools.product(range(3), repeat=n):
        Pm = torch.ones(1, 1, dtype=cplx)
        for m in mus: Pm = torch.kron(Pm, PAULI[m])
        fullP.append(Pm)
    fullP = torch.stack(fullP)
    def q_qubit(X): return (torch.abs(torch.einsum('pij,ji->p', fullP, X)) ** 2).sum() / 2 ** n
    def member(p):
        m = 2 ** n + 1
        A = torch.complex(p[:2 * m], p[2 * m:4 * m]).reshape(m, 2); B = torch.complex(p[4 * m:6 * m], p[6 * m:]).reshape(m, 2)
        Af = torch.zeros(D, 2, dtype=cplx); Bf = torch.zeros(D, 2, dtype=cplx)
        Af[Kidx] = A[:-1]; Af[eidx] = A[-1]; Bf[Kidx] = B[:-1]; Bf[eidx] = B[-1]
        C = Af @ Bf.conj().T; return C / torch.linalg.norm(C)
    return D, Phi, q, Kidx, eidx, q_qubit, member
def lbfgs(f, p, iters=2000):
    p = p.clone().requires_grad_(True)
    opt = torch.optim.LBFGS([p], lr=1, max_iter=iters, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = f(p); v.backward(); return v
    opt.step(closure); return p.detach()
for n in range(2, NMAX + 1):
    D, Phi, q, Kidx, eidx, q_qubit, member = make_tools(n); m = 2 ** n + 1
    err = 0.0
    for k in range(20):
        C = member(torch.randn(8 * m))
        CK = C[Kidx][:, Kidx]; kap = C[Kidx, eidx]; kapp = C[eidx, Kidx]; g = C[eidx, eidx]
        rhs = q_qubit(CK) + torch.linalg.norm(kap) ** 2 + torch.linalg.norm(kapp) ** 2 + abs(g) ** 2 / 2 ** n \
              + 2 * (-0.5) ** n * (g.conj() * torch.trace(CK)).real
        err = max(err, abs(q(C) - rhs).item())
    def gap(p):
        X = torch.complex(p[:2 ** n * 2], p[2 ** n * 2:2 ** n * 4]).reshape(2 ** n, 2) @ \
            torch.complex(p[2 ** n * 4:2 ** n * 6], p[2 ** n * 6:]).reshape(2 ** n, 2).conj().T
        X = X / torch.linalg.norm(X); s = torch.linalg.svdvals(X)
        return 2 ** n * q_qubit(X) - torch.relu(torch.abs(torch.trace(X)) - 2 ** n * s[1]) ** 2
    gmin = min(gap(lbfgs(gap, torch.randn(2 ** n * 8))).item() for _ in range(4))
    rmin = min(q(member(lbfgs(lambda p: q(member(p)), torch.randn(8 * m)))).item() for _ in range(4))
    print(f"n = {n}: identity max error {err:.1e} (20 random members); qubit lemma min gap {gmin:+.1e} (4 adversarial runs); "
          f"min q/||C||^2 over the class {rmin:+.1e} (4 adversarial runs)", flush=True)
# (4) the two-copy critical point
C = np.load(__file__.replace('check_one_excitation.py', 'critical_n2.npy')); C = C / np.linalg.norm(C)
U, S, Vh = np.linalg.svd(C); A = U[:, :2]
T = (A @ A.conj().T).reshape(3, 3, 3, 3)
G = np.kron(np.linalg.eigh(np.einsum('ajbj->ab', T))[1], np.linalg.eigh(np.einsum('jajb->ab', T))[1])   # one local unitary
Cn = G.conj().T @ C @ G
odd = [i for i in range(9) if ((i // 3) == 1) != ((i % 3) == 1)]       # strings with exactly one middle index
even_K = [i for i in range(9) if (i // 3) != 1 and (i % 3) != 1]       # K (x) K with K = span of eigenvectors 0, 2
R2 = float(make_tools(2)[2](torch.tensor(C)))
print(f"two-copy critical point: R = {R2:.6f} "
      f"(exact (9 - sqrt 41)/20 = {(9 - np.sqrt(41)) / 20:.6f}); weight of Cn outside (K (x) K + C|11>) on either side: "
      f"{np.linalg.norm(Cn[odd]) ** 2 + np.linalg.norm(Cn[:, odd]) ** 2:.1e}")
