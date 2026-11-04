"""Brute force on the class of the theorem: rank-two X with range in H_1 (single excitations over |0^n>) and co-range in
H_0 + H_1 + H_2. Random and adversarial minimisation of q_n(X)/||X||^2 by direct partial traces, n = 2..5.
Also, for the open sectors (co-range in H_3), the same search at n = 4, 5 (numerical evidence only)."""
import numpy as np
from scipy.optimize import minimize
from common import q_direct, strings, weight, single_excitations, embed

rng = np.random.default_rng(5)


def search(n, allowed, starts, rand):
    rows = single_excitations(n)
    cols = [st for st in strings(n) if weight(st) in allowed]
    R, C = len(rows), len(cols)
    def build(x):
        U = (x[:2 * R] + 1j * x[2 * R:4 * R]).reshape(R, 2)
        V = (x[4 * R:4 * R + 2 * C] + 1j * x[4 * R + 2 * C:]).reshape(C, 2)
        return embed(U @ V.conj().T, n, rows, cols)
    obj = lambda x: (lambda X: q_direct(X, n) / np.linalg.norm(X) ** 2)(build(x))
    npar = 4 * R + 4 * C
    rmin = min(obj(rng.normal(size=npar)) for _ in range(rand))
    amin = min(minimize(obj, rng.normal(size=npar), method='L-BFGS-B', options={'maxiter': 500}).fun
               for _ in range(starts))
    return rmin, amin


overall = np.inf
for n in range(2, 6):
    rmin, amin = search(n, {0, 1, 2}, 6 if n <= 4 else 3, 300)
    overall = min(overall, amin, rmin)
    print(f"theorem class n={n}: random min q/||X||^2 = {rmin:+.3e}, adversarial min = {amin:+.3e}", flush=True)
for n in (4, 5):
    rmin, amin = search(n, {3}, 4 if n == 4 else 2, 200)
    print(f"open sector c=3, n={n} (evidence only): random min = {rmin:+.3e}, adversarial min = {amin:+.3e}", flush=True)
print(f"RESULT: smallest value on the theorem class {overall:+.3e}")
assert overall > -1e-6
print("PASS")
