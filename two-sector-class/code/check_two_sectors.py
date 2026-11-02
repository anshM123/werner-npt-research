"""Two orthogonal product sectors with random orthogonal subspaces A_s, B_s (dims (2,1), (1,2) or (1,1)) in random
bases.  Checks the identity
  2^n q(C) = 2^a Q_a(X) + 2^b Q_b(Y) + 2^n (|C_AB|^2 + |C_BA|^2) + 2 (-1)^n Re(Tr X conj Tr Y)
(a, b = numbers of copies where A_s resp. B_s is 2-dim; Q = (1/2^k) sum over full-weight Pauli strings |Tr P X|^2 on the
2-dim factors, 1-dim factors trivial) and minimises q/||C||^2 over the class.  usage: python check_two_sectors.py n trials seed"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, TRIALS, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); D = 3 ** n; torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
cplx = torch.complex128
w9 = torch.eye(3, dtype=cplx).reshape(9) / np.sqrt(3)
def Phi(C):
    X = C.reshape((3,) * (2 * n)).permute([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n)
    for s in range(n):
        tr = torch.tensordot(w9.conj(), X, dims=([0], [s])); X = X - 1.5 * torch.moveaxis(torch.tensordot(w9, tr, dims=0), 0, s)
    return X.reshape([3, 3] * n).permute([2 * s for s in range(n)] + [2 * s + 1 for s in range(n)]).reshape(D, D)
def q(C): return torch.vdot(C.reshape(-1), Phi(C).reshape(-1)).real
PA = [torch.tensor(m, dtype=cplx) for m in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])]
def Qk(X, k):
    if k == 0: return torch.abs(X.reshape(())) ** 2 if X.numel() == 1 else None
    tot = 0
    for mus in itertools.product(range(3), repeat=k):
        Pm = torch.ones(1, 1, dtype=cplx)
        for m in mus: Pm = torch.kron(Pm, PA[m])
        tot = tot + torch.abs(torch.trace(Pm @ X)) ** 2
    return tot / 2 ** k
worst_id, worst_min = 0.0, 1e9
for trial in range(TRIALS):
    dims = [[(2, 1), (1, 2), (1, 1)][rng.integers(3)] for _ in range(n)]
    Us = []
    for s in range(n):
        Z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); Q_, _ = np.linalg.qr(Z); Us.append(torch.tensor(Q_))
    Abas = [Us[s][:, :dims[s][0]] for s in range(n)]; Bbas = [Us[s][:, dims[s][0]:dims[s][0] + dims[s][1]] for s in range(n)]
    def kron_all(mats):
        out = torch.ones(1, 1, dtype=cplx)
        for M in mats: out = torch.kron(out, M)
        return out
    EA = kron_all(Abas); EB = kron_all(Bbas); E = torch.cat([EA, EB], 1); mA, mB = EA.shape[1], EB.shape[1]; m = mA + mB
    a = sum(d[0] == 2 for d in dims); b = sum(d[1] == 2 for d in dims)
    def member(p):
        F = torch.complex(p[:2 * m], p[2 * m:4 * m]).reshape(m, 2); G = torch.complex(p[4 * m:6 * m], p[6 * m:]).reshape(m, 2)
        Cs = F @ G.conj().T; Cs = Cs / torch.linalg.norm(Cs); return E @ Cs @ E.conj().T, Cs
    C, Cs = member(torch.randn(8 * m))
    X, Y, P, Qm = Cs[:mA, :mA], Cs[mA:, mA:], Cs[:mA, mA:], Cs[mA:, :mA]
    rhs = 2 ** a * Qk(X, a) + 2 ** b * Qk(Y, b) + 2 ** n * (torch.linalg.norm(P) ** 2 + torch.linalg.norm(Qm) ** 2) \
          + 2 * (-1) ** n * (torch.trace(X) * torch.trace(Y).conj()).real
    worst_id = max(worst_id, abs(2 ** n * q(C) - rhs).item())
    best = 1e9
    for st in range(3):
        p = torch.randn(8 * m, requires_grad=True)
        opt = torch.optim.LBFGS([p], lr=1, max_iter=2000, tolerance_grad=1e-14, tolerance_change=1e-18, history_size=50, line_search_fn='strong_wolfe')
        def closure():
            opt.zero_grad(); v = q(member(p)[0]); v.backward(); return v
        opt.step(closure); best = min(best, q(member(p.detach())[0]).item())
    worst_min = min(worst_min, best)
    print(f"trial {trial}: dims {dims}, identity error {abs(2 ** n * q(C) - rhs).item():.1e}, min q/||C||^2 {best:+.2e}", flush=True)
print(f"n = {n}: max identity error {worst_id:.1e}; min over all trials {worst_min:+.2e}")
