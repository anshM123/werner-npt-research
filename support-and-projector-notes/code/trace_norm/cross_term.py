"""Section 9: the cross term between orthogonal pure states. For orthonormal psi1, psi2 in (C^3)^n, X = psi1 psi2^*, is |q_n(2, X)| <= 1 ?
q_n(2,X) = sum_M (-2)^|M| ||Tr_M X||^2.  This equals (-1)^n 2^n <rho1, Phi^{(x)n} rho2> (Phi at t = 1/2).
usage: python cross_term.py n runs seed sign   (sign=+1 minimises (-1)^n q, i.e. tests the PSD bound; -1 maximises)"""
import sys, itertools
import numpy as np
import torch
torch.set_default_dtype(torch.float64)
n, runs, seed, sg = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
d = 3 ** n
subsets = [K for r in range(n + 1) for K in itertools.combinations(range(n), r)]

def q2(X):
    T = X.reshape([3] * n + [3] * n)
    tot = 0.
    for K in subsets:
        Y, rem = T, list(range(n))
        for s in sorted(K, reverse=True):
            i = rem.index(s)
            Y = torch.diagonal(Y, dim1=i, dim2=len(rem) + i).sum(-1)
            rem.pop(i)
        tot = tot + (-2.0) ** len(K) * (Y.real ** 2 + Y.imag ** 2).sum()
    return tot

def val(params):
    Z = torch.complex(params[0], params[1])          # d x 2
    Q, _ = torch.linalg.qr(Z)                        # orthonormal pair
    X = torch.outer(Q[:, 0], Q[:, 1].conj())
    return sg * (-1) ** n * q2(X)

rng = np.random.default_rng(seed)
out = []
for r in range(runs):
    params = [torch.tensor(rng.standard_normal((d, 2)), requires_grad=True) for _ in range(2)]
    opt = torch.optim.LBFGS(params, lr=1, max_iter=3000, tolerance_grad=1e-13, tolerance_change=1e-16,
                            history_size=50, line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad(); v = val(params); v.backward(); return v
    for _ in range(3):
        opt.step(closure)
    out.append(val(params).item())
print(f"n={n} sign={sg}: min of {'(-1)^n q' if sg > 0 else '-(-1)^n q'} = {min(out):.6f}  (B_n needs >= -1)")
