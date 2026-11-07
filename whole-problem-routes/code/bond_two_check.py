"""Floating-point illustration of Proposition 7.17 (the proposition itself is proved by hand).

For open matrix product states x of bond dimension two the proposition gives lambda_min(P_a^{(x)n}(x x^*)) >= (a - 1/2)^n
|x|^2.  This script computes lambda_min / (a - 1/2)^n at a = 3/4 for
  - the open AKLT chain (Cartesian spin-1 basis, A^k = sigma_k / sqrt 3, random boundary vectors), whose one-copy
    reduced states have full rank away from the ends, so Corollary 7.14 does not cover it;
  - random complex open matrix product states of bond dimension two;
  - for comparison, random ones of bond dimension three (not covered by the proposition).
Usage: python bond_two_check.py n samples seed"""
import numpy as np, sys
n, samples, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(seed)
a = 0.75
S = np.zeros((9, 9))                       # P_a as a matrix on vec(X), X[i, j] -> i*3 + j
for i in range(3):
    for j in range(3):
        out = np.zeros((3, 3))
        if i == j:
            out += a * np.eye(3)
        out[i, j] -= 0.25
        out[j, i] -= 0.25
        S[:, i*3 + j] = out.ravel()

def lam_min(x):
    X = np.outer(x, x.conj()).reshape([3] * (2 * n))
    for s in range(n):
        X = np.moveaxis(X, [s, n + s], [0, 1])
        sh = X.shape
        X = (S @ X.reshape(9, -1)).reshape(sh)
        X = np.moveaxis(X, [0, 1], [s, n + s])
    M = X.reshape(3**n, 3**n)
    return np.linalg.eigvalsh((M + M.conj().T) / 2)[0]

def chain(tensors, left, right):
    v = left.reshape(1, -1)                 # shape (1, D)
    for A in tensors:                       # A: (3, D, D)
        v = np.einsum('...a,kab->...kb', v, A)
    v = np.einsum('...b,b->...', v, right).reshape(-1)
    return v / np.linalg.norm(v)

def marg_ranks(x):
    T = x.reshape([3] * n); out = []
    for s in range(n):
        M = np.moveaxis(T, s, 0).reshape(3, -1)
        out.append(int(np.sum(np.linalg.svd(M, compute_uv=False) > 1e-9)))
    return out

sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
A_aklt = np.array(sig) / np.sqrt(3)
prod = (a - 0.5)**n
print(f"n = {n}, a = 3/4, product value (a - 1/2)^n = {prod:.6g}")
vals = []
for t in range(samples):
    l = rng.standard_normal(2) + 1j * rng.standard_normal(2); r = rng.standard_normal(2) + 1j * rng.standard_normal(2)
    x = chain([A_aklt] * n, l, r)
    vals.append(lam_min(x) / prod)
print(f"open AKLT chain, random boundary vectors: min ratio {min(vals):.4f}; one-copy ranks {marg_ranks(x)}")
for D in (2, 3):
    vals = []
    for t in range(samples):
        T = [rng.standard_normal((3, D, D)) + 1j * rng.standard_normal((3, D, D)) for _ in range(n)]
        l = rng.standard_normal(D) + 1j * rng.standard_normal(D); r = rng.standard_normal(D) + 1j * rng.standard_normal(D)
        vals.append(lam_min(chain(T, l, r)) / prod)
    print(f"random open MPS, bond dimension {D}: min ratio {min(vals):.4f} over {samples}")
print("ratios >= 1 expected for the AKLT chain and bond dimension two (Proposition 7.17)")
