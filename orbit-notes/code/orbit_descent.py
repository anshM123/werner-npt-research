"""Descent of the Werner form along local-filter (squeeze) orbits.

For an operator C of rank at most two on (C^3)^{(x)n} and positive definite G_s = expm(H_s) at each copy, minimise
    R_t(G C G),   G = G_1 (x) ... (x) G_n,   R_t(X) = sum_J (-t)^|J| ||Tr_J X||^2 / ||X||^2,
over the Hermitian H_s by L-BFGS.  The run reports where the descent stops: at the boundary of the orbit (the operator
becomes thin at some copy, i.e. the third eigenvalue of a local support goes to 0, or its second singular value goes to
0), or at an interior point.  An interior stop would be a local minimum of R_t on the orbit.

usage: python orbit_descent.py n t starts seed ensemble
  ensemble: random | projector | psd | hermitian | nearzero"""
import sys, numpy as np, torch
from scipy.optimize import minimize
torch.set_default_dtype(torch.float64)
n, t, STARTS, SEED, ens = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
D = 3 ** n
torch.manual_seed(SEED); rng = np.random.default_rng(SEED)
e = torch.eye(3, dtype=torch.complex128).reshape(9) / np.sqrt(3)

def pairs(C):
    T = C.reshape((3,) * (2 * n)); perm = [i for s in range(n) for i in (s, n + s)]
    return T.permute(perm).reshape((9,) * n)

def R(C):
    """R_t(C): apply I - 3t |e><e| (e = vec(I)/sqrt3) at every copy of vec C."""
    X = pairs(C); Y = X
    for s in range(n):
        tr = torch.tensordot(e.conj(), Y, dims=([0], [s]))
        Y = Y - 3 * t * torch.moveaxis(torch.tensordot(e, tr, dims=0), 0, s)
    return torch.vdot(X.reshape(-1), Y.reshape(-1)).real / torch.vdot(X.reshape(-1), X.reshape(-1)).real

def herm(x):
    M = torch.zeros(3, 3, dtype=torch.complex128) + torch.diag(x[:3].to(torch.complex128))
    iu = torch.triu_indices(3, 3, 1)
    M[iu[0], iu[1]] = torch.complex(x[3:6], x[6:9]); M[iu[1], iu[0]] = torch.complex(x[3:6], -x[6:9])
    return M

def squeeze(C, x):
    G = torch.ones(1, 1, dtype=torch.complex128)
    for s in range(n):
        G = torch.kron(G, torch.linalg.matrix_exp(herm(x[9 * s:9 * s + 9])))
    return G @ C @ G

def thickness(C):
    """(smallest third eigenvalue of the local supports of range and co-range over all copies, sigma_2/sigma_1)"""
    U, S, Vh = torch.linalg.svd(C); out = []
    for M in (U[:, :2], Vh[:2].conj().T):
        P = (M @ M.conj().T).reshape((3,) * (2 * n)); L = 'abcdefghijklmnopqrstuvwxyz'
        for s in range(n):
            ins = [L[i] for i in range(2 * n)]
            for a in range(n):
                if a != s: ins[n + a] = ins[a]
            out.append(torch.linalg.eigvalsh(torch.einsum(''.join(ins) + '->' + L[s] + ins[n + s], P))[0].item())
    return min(out), (S[1] / S[0]).item()

def start():
    cg = lambda *sh: torch.randn(*sh, dtype=torch.complex128)
    if ens == 'random':
        return cg(D, 2) @ cg(D, 2).conj().T
    if ens == 'projector':
        Q, _ = torch.linalg.qr(cg(D, 2)); return Q @ Q.conj().T
    if ens == 'psd':
        A = cg(D, 2); return A @ A.conj().T
    if ens == 'hermitian':
        Q, _ = torch.linalg.qr(cg(D, 2))
        return Q[:, :1] @ Q[:, :1].conj().T + float(rng.uniform(-1, 1)) * Q[:, 1:] @ Q[:, 1:].conj().T
    if ens == 'nearzero':       # thin zero P_V (x) |x><y| plus a thick perturbation, projected back to rank two
        x = cg(3 ** (n - 1), 1); y = cg(3 ** (n - 1), 1); U, _ = torch.linalg.qr(cg(3, 3))
        P = U @ torch.diag(torch.tensor([1, 1, 0], dtype=torch.complex128)) @ U.conj().T
        C0 = torch.kron(P, x @ y.conj().T); C0 = C0 / torch.linalg.norm(C0)
        E = cg(D, 2) @ cg(D, 2).conj().T; E = E / torch.linalg.norm(E)
        Uc, Sv, Vh = torch.linalg.svd(C0 + float(rng.choice([1e-3, 1e-2, 0.1])) * E)
        return (Uc[:, :2] * Sv[:2]) @ Vh[:2]
    raise SystemExit('unknown ensemble')

print(f"n = {n}, t = {t}, ensemble {ens}; product value (1-2t)(1-t)^(n-1) = {(1 - 2 * t) * (1 - t) ** (n - 1):+.6f}, "
      f"rank-one value (1-t)^n = {(1 - t) ** n:+.6f}")
interior = 0
for k in range(STARTS):
    C = start(); C = C / torch.linalg.norm(C)
    def f(xn):
        x = torch.tensor(xn, requires_grad=True); r = R(squeeze(C, x)); r.backward()
        return r.item(), x.grad.numpy().copy()
    r0 = R(C).item(); th0, _ = thickness(C)
    res = minimize(f, np.zeros(9 * n), jac=True, method='L-BFGS-B', options={'maxiter': 4000, 'gtol': 1e-11})
    th, sr = thickness(squeeze(C, torch.tensor(res.x)))
    inside = th > 1e-3 and sr > 1e-3; interior += inside
    print(f"  start {k}: R {r0:+.6f} (thickness {th0:.2e}) -> {res.fun:+.6f}; end thickness {th:.1e}, s2/s1 {sr:.1e}, "
          f"max|H| {np.abs(res.x).max():.1f}" + ("   INTERIOR STOP" if inside else ""), flush=True)
print(f"interior stops: {interior} of {STARTS}")
