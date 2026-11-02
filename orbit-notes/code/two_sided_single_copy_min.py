"""Check data/two_sided_single_copy_min_n3.npy (n = 3, t = 1/2): a thick rank-two C that is a STRICT local minimum of R
over its two-sided single-copy orbit at copy 0,  C -> (e^X (x) 1 (x) 1) C (e^Y (x) 1 (x) 1)^T,  X, Y in gl(3).
Prints R, the thickness at every copy, the gradient over (X, Y) at each copy, the Hessian over (X, Y) at copy 0 on the
complement of the directions along which R is constant (unitary conjugations and scalings; their span has rank 12,
found with a rank-revealing SVD), and R along 300 random directions of that complement.
usage: python two_sided_single_copy_min.py   (run from orbit-notes/)"""
import numpy as np
import torch

torch.set_default_dtype(torch.float64)
n, t = 3, 0.5
D = 3 ** n
e9 = torch.eye(3, dtype=torch.complex128).reshape(9) / np.sqrt(3)


def R_t(C):
    X = C.reshape((3,) * (2 * n)).permute([i for s in range(n) for i in (s, n + s)]).reshape((9,) * n); Y = X
    for s in range(n):
        tr = torch.tensordot(e9.conj(), Y, dims=([0], [s]))
        Y = Y - 3 * t * torch.moveaxis(torch.tensordot(e9, tr, dims=0), 0, s)
    return torch.vdot(X.reshape(-1), Y.reshape(-1)).real / torch.vdot(X.reshape(-1), X.reshape(-1)).real


def emb(M, s):
    out = torch.ones(1, 1, dtype=torch.complex128)
    for k in range(n):
        out = torch.kron(out, M if k == s else torch.eye(3, dtype=torch.complex128))
    return out


def act(C, z, s):
    X = torch.complex(z[:9], z[9:18]).reshape(3, 3); Y = torch.complex(z[18:27], z[27:]).reshape(3, 3)
    return emb(torch.linalg.matrix_exp(X), s) @ C @ emb(torch.linalg.matrix_exp(Y), s).T


def invariant_dirs():
    dirs = []
    for M in (torch.eye(3, dtype=torch.complex128), 1j * torch.eye(3, dtype=torch.complex128)):
        dirs.append(torch.cat([M.real.reshape(-1), M.imag.reshape(-1), torch.zeros(18)]))
        dirs.append(torch.cat([torch.zeros(18), M.real.reshape(-1), M.imag.reshape(-1)]))
    for a in range(3):
        for b in range(a, 3):
            for ph in ((1.0,) if a == b else (1.0, 1j)):
                K = torch.zeros(3, 3, dtype=torch.complex128); K[a, b] = ph; K[b, a] = np.conj(ph)
                X = 1j * K; Y = -1j * K.T
                dirs.append(torch.cat([X.real.reshape(-1), X.imag.reshape(-1), Y.real.reshape(-1), Y.imag.reshape(-1)]))
    N = torch.stack(dirs, 1)
    U, S, _ = torch.linalg.svd(N, full_matrices=False)
    return U[:, S > 1e-10 * S[0]]


def thickness_each(C):
    out = []
    for M in (C, C.conj().T):
        U, S, Vh = torch.linalg.svd(M); P = (U[:, :2] @ U[:, :2].conj().T).reshape((3,) * (2 * n))
        for s in range(n):
            T = torch.movedim(P, (s, n + s), (0, 1)).reshape(3, 3, 3 ** (n - 1), 3 ** (n - 1))
            out.append(torch.linalg.eigvalsh(torch.einsum('abxx->ab', T))[0].item())
    return out


C = torch.tensor(np.load('data/two_sided_single_copy_min_n3.npy')); C = C / torch.linalg.norm(C)
r0 = R_t(C).item()
print(f"R = {r0:+.8f}, singular values {np.round(torch.linalg.svdvals(C)[:3].numpy(), 4)}")
print(f"thickness (range at copies 0,1,2; co-range at copies 0,1,2): {np.round(thickness_each(C), 4)}")
for s in range(n):
    z = torch.zeros(36, requires_grad=True)
    g = torch.autograd.grad(R_t(act(C, z, s)), z)[0]
    print(f"copy {s}: |gradient over (X, Y)| {torch.linalg.norm(g).item():.1e}")
f = lambda zz: R_t(act(C, zz, 0))
H = torch.autograd.functional.hessian(f, torch.zeros(36)); H = (H + H.T) / 2
Q = invariant_dirs()
print(f"directions along which R is constant: rank {Q.shape[1]}; largest change of R along them at step 1e-2: "
      f"{max(abs(f(1e-2 * Q[:, k]).item() - r0) for k in range(Q.shape[1])):.1e}")
P = torch.eye(36) - Q @ Q.T; w, U = torch.linalg.eigh(P); B = U[:, w > 0.5]
ev = torch.linalg.eigvalsh(B.T @ H @ B)
print(f"Hessian at copy 0 on the {B.shape[1]}-dim complement: smallest eigenvalues {np.round(ev[:4].numpy(), 5)}")
rng = np.random.default_rng(0)
for eps in (1e-3, 1e-2, 5e-2):
    vals = []
    for k in range(300):
        d = B @ torch.tensor(rng.standard_normal(B.shape[1])); d = d / torch.linalg.norm(d)
        vals.append(f(eps * d).item() - r0)
    print(f"step {eps:g}: smallest R(step) - R over 300 random directions of the complement {min(vals):+.2e}")
