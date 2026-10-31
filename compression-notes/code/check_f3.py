"""Section 11, F3: is every 3 x 3 local filtering of rho^{(x)m} in K_1(3) = cone{products, (a (x) b) rho (a (x) b)^*}?

Targets: start from the single-copy filtering (all copies but the first contracted with |0>, which gives rho itself,
partial-transpose minimum -1/15) and move the filters by delta toward random m-copy filters.  For delta > 0 the target
is a genuinely m-copy filtering; for small delta it is NPT.  Each target is fitted by J filterings of rho plus P
products (relative Frobenius residual).  Controls: a random element of K_1(3), and the antisymmetric projector (it is
1-distillable, so it lies outside).  usage: python check_f3.py m   (m = 3 or 4; needs PyTorch, about 15 minutes)
"""
import itertools
import sys

import numpy as np
import torch

torch.set_default_dtype(torch.float64)
m = int(sys.argv[1]) if len(sys.argv) > 1 else 3
J, P, RS = 6, 12, 4
D = 3
subsets = [S for r in range(m + 1) for S in itertools.combinations(range(m), r)]


def filt(A, B):
    """(A (x) B) rho^{(x)m} (A (x) B)^*, using rho^{(x)m} = sum_S (-1/2)^{|S|} F_S."""
    At, Bt = A.reshape([D] + [3] * m), B.reshape([D] + [3] * m)
    Y = 0
    for S in subsets:
        Sc = [s for s in range(m) if s not in S]
        pa = [0] + [1 + s for s in S] + [1 + s for s in Sc]
        Ap = At.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc))
        Bp = Bt.permute(pa).reshape(D, 3 ** len(S), 3 ** len(Sc))
        M1 = torch.einsum('axc,eyc->aexy', Ap, Ap.conj())
        M2 = torch.einsum('byc,fxc->bfyx', Bp, Bp.conj())
        Y = Y + (-0.5) ** len(S) * torch.einsum('aexy,bfyx->abef', M1, M2).reshape(D * D, D * D)
    return Y


F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
R1 = torch.tensor(np.eye(9) - F / 2, dtype=torch.complex128)


def pt(Y):
    return Y.reshape(3, 3, 3, 3).permute(0, 3, 2, 1).reshape(9, 9)


def fit(Y):
    Yt = torch.tensor(Y / np.trace(Y).real, dtype=torch.complex128)
    best = np.inf
    for _ in range(RS):
        p = [(0.4 * torch.randn(*s, dtype=torch.complex128)).requires_grad_()
             for s in ((J, 3, 3), (J, 3, 3), (P, 3), (P, 3))]

        def model():
            V = torch.einsum('ka,kb->kab', p[2], p[3]).reshape(P, 9)
            Kt = torch.einsum('jab,jcd->jacbd', p[0], p[1]).reshape(J, 9, 9)
            return V.T @ V.conj() + (Kt @ R1 @ Kt.conj().transpose(1, 2)).sum(0)
        opt = torch.optim.LBFGS(p, lr=1, max_iter=6000, history_size=50, line_search_fn='strong_wolfe',
                                tolerance_grad=1e-15, tolerance_change=1e-18)

        def closure():
            opt.zero_grad()
            loss = torch.linalg.norm(model() - Yt) ** 2
            loss.backward()
            return loss
        opt.step(closure)
        with torch.no_grad():
            best = min(best, float(torch.linalg.norm(model() - Yt) / torch.linalg.norm(Yt)))
    return best


torch.manual_seed(0)
g = np.random.default_rng(1)
a = g.standard_normal((2, 3, 3)) + 1j * g.standard_normal((2, 3, 3))
b = g.standard_normal((2, 3, 3)) + 1j * g.standard_normal((2, 3, 3))
Yc = sum(np.kron(a[i], b[i]) @ R1.numpy() @ np.kron(a[i], b[i]).conj().T for i in range(2))
for _ in range(5):
    x, y = g.standard_normal(3) + 1j * g.standard_normal(3), g.standard_normal(3) + 1j * g.standard_normal(3)
    v = np.kron(x, y)
    Yc = Yc + np.outer(v, v.conj())
print(f"control in K_1(3): fit {fit(Yc):.2e}", flush=True)
print(f"control antisymmetric projector (outside): fit {fit((np.eye(9) - F) / 2):.2e}", flush=True)

u = torch.zeros(3 ** (m - 1), dtype=torch.complex128)
u[0] = 1
A0 = torch.einsum('ax,y->axy', torch.eye(3, dtype=torch.complex128), u).reshape(3, 3 ** m)
g2 = torch.Generator().manual_seed(5)
RA = torch.randn(3, 3 ** m, dtype=torch.complex128, generator=g2)
RB = torch.randn(3, 3 ** m, dtype=torch.complex128, generator=g2)
for delta in [0.0, 0.1, 0.2, 0.3]:
    A = A0 / torch.linalg.norm(A0) + delta * RA / torch.linalg.norm(RA)
    B = A0 / torch.linalg.norm(A0) + delta * RB / torch.linalg.norm(RB)
    with torch.no_grad():
        Y = filt(A, B)
        Y = (Y + Y.conj().T) / 2
        Y = Y / torch.trace(Y).real
        f = torch.linalg.eigvalsh(pt(Y))[0].item()
    print(f"m = {m}, delta {delta}: PT-min {f:+.4f}, K_1(3) fit {fit(Y.numpy()):.2e}", flush=True)
