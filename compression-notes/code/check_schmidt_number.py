"""Remark (Schmidt number of two and three copies), Section 9.

Checks the identities
    rho^{(x)2}/Tr = 0.24 W(-1)(x)W(0) + 0.24 W(0)(x)W(-1) + 0.04 W(1)(x)W(1) + 0.48 W(0)(x)W(0),
    rho^{(x)3}/Tr = 0.288 [W(-1)W(0)W(0) + perms] + 0.024 [W(-1)W(1)W(1) + perms] + 0.064 W(1)W(1)W(1),
with W(s) the normalised Werner state on C^3 (x) C^3 with <F> = s (rho = I - F/2 is W(-1/5) up to normalisation),
and that antisymmetric vectors of C^3 (x) C^3 have Schmidt rank 2 (so W(-1) = P_a/3 has Schmidt number 2).
Then fits rho (x) rho by sums of Schmidt-rank-2 pure terms, with controls (needs PyTorch; about a minute).
"""
import itertools
import numpy as np

F = np.zeros((9, 9))
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
I9 = np.eye(9)
Ps, Pa = (I9 + F) / 2, (I9 - F) / 2
rho = I9 - F / 2


def W(s):
    return (1 + s) / 2 * Ps / 6 + (1 - s) / 2 * Pa / 3


def copies(ops):
    """Tensor product of per-copy operators, reordered to Alice = (copies), Bob = (copies)."""
    m = len(ops)
    R = np.eye(1)
    for o in ops:
        R = np.kron(R, o)
    perm = [2 * i for i in range(m)] + [2 * i + 1 for i in range(m)]
    return R.reshape([3] * (4 * m)).transpose(perm + [2 * m + q for q in perm]).reshape(9 ** m, 9 ** m)


print("Tr(W(s) F) for s = -1, -0.2, 0, 1:", [round(float(np.trace(W(s) @ F)), 12) for s in (-1, -0.2, 0, 1)])
print("rho/Tr rho = W(-1/5):", np.allclose(rho / np.trace(rho), W(-0.2)))

target2 = copies([rho, rho]) / np.trace(rho) ** 2
mix2 = (0.24 * copies([W(-1), W(0)]) + 0.24 * copies([W(0), W(-1)]) + 0.04 * copies([W(1), W(1)])
        + 0.48 * copies([W(0), W(0)]))
print(f"two copies: max |difference| = {np.abs(target2 - mix2).max():.1e}")

target3 = copies([rho, rho, rho]) / np.trace(rho) ** 3
mix3 = 0.064 * copies([W(1)] * 3)
for p in set(itertools.permutations([-1, 0, 0])):
    mix3 = mix3 + 0.288 * copies([W(s) for s in p])
for p in set(itertools.permutations([-1, 1, 1])):
    mix3 = mix3 + 0.024 * copies([W(s) for s in p])
print(f"three copies: max |difference| = {np.abs(target3 - mix3).max():.1e}")

g = np.random.default_rng(0)
ranks = []
for _ in range(1000):
    M = g.standard_normal((3, 3)) + 1j * g.standard_normal((3, 3))
    ranks.append(np.linalg.matrix_rank(M - M.T, 1e-10))
print("Schmidt ranks of 1000 random antisymmetric vectors:", sorted(set(ranks)))

try:
    import torch
except ImportError:
    torch = None
if torch is not None:
    torch.set_default_dtype(torch.float64)

    def fit(X, P=300, IT=4000, RS=2):
        Xt = torch.tensor(X / np.trace(X).real, dtype=torch.complex128)
        best = np.inf
        for _ in range(RS):
            p = [(0.1 * torch.randn(P, 9, dtype=torch.complex128)).requires_grad_() for _ in range(4)]

            def model():
                V = (torch.einsum('ka,kb->kab', p[0], p[1]) + torch.einsum('ka,kb->kab', p[2], p[3])).reshape(P, 81)
                return V.T @ V.conj()
            opt = torch.optim.LBFGS(p, lr=1, max_iter=IT, history_size=50, line_search_fn='strong_wolfe',
                                    tolerance_grad=1e-15, tolerance_change=1e-18)

            def closure():
                opt.zero_grad()
                loss = torch.linalg.norm(model() - Xt) ** 2
                loss.backward()
                return loss
            opt.step(closure)
            with torch.no_grad():
                best = min(best, float(torch.linalg.norm(model() - Xt) / torch.linalg.norm(Xt)))
        return best
    phi = np.eye(9).reshape(81) / 3
    Xe = 0.5 * np.outer(phi, phi) + 0.5 * np.eye(81) / 81
    print(f"fit by Schmidt-rank-2 terms: control with Schmidt number > 2 (fidelity 1/2 with Phi_9): {fit(Xe):.2e}")
    print(f"fit by Schmidt-rank-2 terms: rho (x) rho: {fit(copies([rho, rho])):.2e}")
