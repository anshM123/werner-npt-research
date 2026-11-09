"""Local minimisation of q_N(1/2, C)/||C||_F^2 over rank-two C = A B^* (complex), started at the star near-zeros of
Section 7 (thick at every copy, q/||C||^2 = 2^{-(k+1)} when k + 1 = m) plus a small perturbation.
usage: python star_descent.py k m eps steps seed   (needs PyTorch)"""
import sys, itertools
import numpy as np
import torch
torch.set_default_dtype(torch.float64)
k, m, eps, steps, seed = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
torch.manual_seed(seed); rng = np.random.default_rng(seed)
n = k + m; N = n + 1; D = 3 ** N; t = 0.5
om = np.exp(2j * np.pi / 3)
X = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)
G = np.zeros((N, N), dtype=np.int64); G[0, 1:] = 1; G[1:, 0] = 1
ph = sum(G[i, j] * X[:, i] * X[:, j] for i in range(N) for j in range(i + 1, N))
g = om ** (ph % 3) / np.sqrt(D)
d = np.array([0] + [0] * k + [1] * m); v = np.zeros(N, dtype=np.int64); v[0] = 1
Zv = lambda u: om ** ((X @ u) % 3)
a1, b1 = Zv(d) * g, g
a2, b2 = Zv((d - v) % 3) * g, Zv(-v % 3) * g
# optimal phase from the closed form: cross term ~ x + w + w^2 y -> choose the phase numerically below
subsets = list(itertools.product([0, 1], repeat=N))
def q_of(A, B):
    # C = A B^*, A, B: D x 2 complex tensors
    tot = 0.0
    TA = A.reshape((3,) * N + (2,)); TB = B.reshape((3,) * N + (2,))
    for K in subsets:
        kept = [i for i in range(N) if not K[i]]; tr = [i for i in range(N) if K[i]]
        dk = 3 ** len(kept); dt = 3 ** len(tr)
        Ma = TA.permute(*kept, *tr, N).reshape(dk, dt, 2)          # a_i as (kept x traced)
        Mb = TB.permute(*kept, *tr, N).reshape(dk, dt, 2)
        # Tr_K C = sum_i Ma_i Mb_i^*  (kept x kept)
        T = torch.einsum('xai,yai->xy', Ma, Mb.conj())
        tot = tot + (-t) ** len(tr) * (T.abs() ** 2).sum()
    return tot
def normC(A, B):
    C = A @ B.conj().T
    return (C.abs() ** 2).sum()
A0 = torch.tensor(np.stack([a1, a2], 1)); B0 = torch.tensor(np.stack([b1, b2], 1))
best = None
for phi in np.linspace(0, 2 * np.pi, 13)[:-1]:
    A = A0.clone(); A[:, 1] = A[:, 1] * np.exp(1j * phi)
    val = (q_of(A, B0) / normC(A, B0)).item()
    if best is None or val < best[0]:
        best = (val, phi)
print(f"star N={N} (k={k}, m={m}): start value q/||C||^2 = {best[0]:.6e} (phase {best[1]:.3f}); 2^-(k+1) = {2.0 ** -(k + 1):.6e}")
A = A0.clone(); A[:, 1] = A[:, 1] * np.exp(1j * best[1])
A = A + eps * torch.randn(D, 2, dtype=torch.complex128); B = B0 + eps * torch.randn(D, 2, dtype=torch.complex128)
ar = torch.view_as_real(A).clone().requires_grad_(True); br = torch.view_as_real(B).clone().requires_grad_(True)
opt = torch.optim.LBFGS([ar, br], lr=0.5, max_iter=steps, tolerance_grad=1e-12, tolerance_change=1e-15, history_size=50,
                        line_search_fn='strong_wolfe')
def closure():
    opt.zero_grad()
    Ac = torch.view_as_complex(ar); Bc = torch.view_as_complex(br)
    val = q_of(Ac, Bc) / normC(Ac, Bc)
    val.backward()
    return val
for it in range(3):
    val = opt.step(closure)
    print(f"  round {it}: q/||C||^2 = {val.item():.6e}")
Ac = torch.view_as_complex(ar).detach(); Bc = torch.view_as_complex(br).detach()
C = (Ac @ Bc.conj().T).numpy()
sv = np.linalg.svd(C, compute_uv=False)[:3]
# local thickness: least eigenvalue of one-copy marginals of range and co-range projectors
U, s, Vh = np.linalg.svd(C); Pr = U[:, :2] @ U[:, :2].conj().T; Pc = Vh[:2].conj().T @ Vh[:2]
def marg_min(Pm):
    T = Pm.reshape((3,) * (2 * N)); out = []
    for c in range(N):
        A_ = T
        cur = N
        for i in sorted([i for i in range(N) if i != c], reverse=True):
            A_ = np.trace(A_, axis1=i, axis2=i + cur); cur -= 1
        out.append(np.linalg.eigvalsh(A_).min())
    return min(out)
print(f"  final singular values {np.round(sv, 5)}, least local eigenvalue range {marg_min(Pr):.4f}, co-range {marg_min(Pc):.4f}")
