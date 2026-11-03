"""Stage 3c: (i) exactness of the determinant-normalised concurrence dual
        C(B) = max(0, sup_{det M = 1} -Tr(W_M B))   vs Wootters' formula;
   (ii) orb24 (2): ||sqrt(S) W sqrt(T)||_* <= sqrt(TrS TrT) + sqrt(Tr(WS) Tr(WT)) for separable S, T, det W-matrix 1;
   (iii) star24 (4): D_x <= ||x||^2 I for D = (I-F_new)(I-F_H);
   (iv) direction  -Tr(W_M B) <= |det M| C(B)  on random B."""
import math, json
import numpy as np
import torch
from core import *

torch.set_num_threads(4)
rng = np.random.default_rng(2026)
out = []


def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)


def Wm(M):
    return np.einsum('kl,ij->kjil', M, M.conj()).reshape(4, 4)


def wit_sup(B, restarts=12, steps=400):
    Bt = torch.tensor(B, dtype=CD)
    best = -1e300
    for r in range(restarts):
        p = torch.randn(8, dtype=torch.float64) * (1 + r % 4)
        p.requires_grad_(True)
        opt = torch.optim.Adam([p], lr=0.05)
        for it in range(steps):
            M = torch.complex(p[:4].reshape(2, 2), p[4:].reshape(2, 2))
            W = torch.einsum('kl,ij->kjil', M, M.conj()).reshape(4, 4)
            val = -(torch.trace(W @ Bt)).real / (torch.linalg.det(M).abs() + 1e-300)
            best = max(best, float(val.detach()))
            opt.zero_grad(); (-val).backward(); opt.step()
    return best


def rand_psd(rank):
    G = rng.normal(size=(4, rank)) + 1j * rng.normal(size=(4, rank))
    return G @ G.conj().T


# (i) dual exactness
maxdev, mindev, n_ent = 0.0, 0.0, 0
for t in range(60):
    rank = 1 + t % 4
    B = rand_psd(rank)
    if t % 5 == 4:   # local filter of a near-separable state
        B = rand_psd(4) + 3 * np.kron(np.diag([1, 0.2]), np.diag([0.3, 1]))
    c = concurrence(B)
    s = wit_sup(B)
    if c > 1e-9:
        n_ent += 1
        maxdev = max(maxdev, (s - c) / np.trace(B).real)
        mindev = min(mindev, (s - c) / np.trace(B).real)
    else:
        maxdev = max(maxdev, max(s, 0) / np.trace(B).real)
log(f"(i) dual: {n_ent} entangled of 60; (sup_witness - Wootters)/TrB in [{mindev:.2e}, {maxdev:.2e}]  (exact dual <=> ~0)")

# (iv) direction: -Tr(W_M B) <= |det M| C(B) for random M
worst = -1e9
for t in range(4000):
    B = rand_psd(1 + t % 4); M = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    worst = max(worst, (-np.trace(Wm(M) @ B).real - abs(np.linalg.det(M)) * concurrence(B)) / (np.trace(B).real * np.linalg.norm(M) ** 2))
log(f"(iv) max(-Tr(W_M B) - |det M| C(B)) normalised = {worst:.2e}  (must be <= 0)")


# (ii) orb24 root-fidelity inequality, adversarial over separable S, T (6 product terms each) and M
def sep_state(p, k=6):
    a = torch.complex(p[:2 * k].reshape(k, 2), p[2 * k:4 * k].reshape(k, 2))
    b = torch.complex(p[4 * k:6 * k].reshape(k, 2), p[6 * k:8 * k].reshape(k, 2))
    v = torch.einsum('ri,rj->rij', a, b).reshape(k, 4)
    return v.T @ v.conj(), p[8 * k:]


def sqrtm_psd(S):
    w, V = torch.linalg.eigh(S)
    return (V * torch.sqrt(torch.clamp(w, min=0.0)).to(CD)) @ V.conj().T


best = -1e9
for r in range(40):
    p = torch.randn(2 * 48 + 8, dtype=torch.float64); p.requires_grad_(True)
    opt = torch.optim.Adam([p], lr=0.02)
    for it in range(600):
        S, q = sep_state(p); T, q = sep_state(q)
        M = torch.complex(q[:4].reshape(2, 2), q[4:8].reshape(2, 2)); M = M / torch.sqrt(torch.linalg.det(M))
        W = torch.einsum('kl,ij->kjil', M, M.conj()).reshape(4, 4)
        lhs = torch.linalg.svdvals(sqrtm_psd(S) @ W @ sqrtm_psd(T)).sum()
        a = torch.trace(W @ S).real; b = torch.trace(W @ T).real
        rhs = torch.sqrt(torch.trace(S).real * torch.trace(T).real) + torch.sqrt(torch.clamp(a * b, min=0.0))
        val = lhs / rhs
        if torch.isfinite(val):
            best = max(best, float(val.detach()))
        else:
            break
        opt.zero_grad(); (-val).backward(); opt.step()
log(f"(ii) orb24 (2): max LHS/RHS over separable S,T = {best:.6f}  (must be <= 1)")

# (iii) D_x <= ||x||^2 I
worst = 0.0
for Dh in (2, 3, 4, 5):
    dims = (2, Dh)
    n = 2 * Dh
    I = torch.eye(n * n, dtype=CD)
    Dop = (I - swap_perm(dims, 0)) @ (I - swap_perm(dims, 1))
    for t in range(300):
        x = rand_c(n).to(CD)
        Dx = torch.einsum('azbw,z,w->ab', Dop.reshape(n, n, n, n), x.conj(), x)   # <x|_2 D |x>_2 on replica 1
        ev = torch.linalg.eigvalsh(0.5 * (Dx + Dx.conj().T))
        worst = max(worst, float(ev.max() / fro2(x)))
        assert float(ev.min()) > -1e-10
log(f"(iii) max eig(D_x)/||x||^2 = {worst:.6f}  (must be <= 1)")
open('stage3c_dual.log', 'w').write('\n'.join(out) + '\n')
