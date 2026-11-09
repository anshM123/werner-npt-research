"""Checks of three all-n statements from an independent audit (Section 7 and OBSTRUCTIONS).
 (A) Complex zeros of the identity-count polynomial approach 0: P_n = |u_n><u_n| + |2^n><2^n|,
     u_n = (sqrt(3^n)|0^n> + sqrt(2^n)|1^n>)/sqrt(3^n + 2^n);  G(x) = 1/2[(2 - B)R^n + B + (2 + B)S^n],
     R = (x+2)/3, S = (x-1)/3, B = 2*6^n/(3^n+2^n)^2 (checked against partial traces, n <= 6), a zero near
     -1/40 + 0.313 i at n = 20, n x_n -> 2 pi i, and the positive endpoint values.
 (B) A real 3x3x3 tensor with negative mixed contraction Tr(rho rho^{T_A}) = -1, its three-copy critical scalar 493,
     and the bound Tr(rho rho^{T_A}) >= (2/d - 1) Tr rho^2 (random tests; equality from quaternion matrices).
 (C) The DDVV construction: local identities and the inequality for n = 1, 2.
usage: python audit89_checks.py"""
import itertools
import numpy as np
from fractions import Fraction as Fr
rng = np.random.default_rng(89)
# ---------------- (A)
def ptrace_norm2(v, n, J):
    T = np.outer(v, v).reshape([3] * (2 * n))
    keep = [s for s in range(n) if s not in J]; L = 'abcdefghij'
    inp = [L[s] for s in range(n)] + [L[s].upper() if s not in J else L[s] for s in range(n)]
    out = [L[s] for s in keep] + [L[s].upper() for s in keep]
    return T, inp, out
def G_direct(n, x):
    p, r = 3 ** n, 2 ** n
    u = np.zeros(3 ** n); u[0] = np.sqrt(p); u[sum(3 ** s for s in range(n))] = np.sqrt(r); u /= np.linalg.norm(u)
    v = np.zeros(3 ** n); v[-1] = 1
    P = (np.outer(u, u) + np.outer(v, v)).reshape([3] * (2 * n)); L = 'abcdefghij'
    tot = 0
    for k in range(n + 1):
        for J in itertools.combinations(range(n), k):
            keep = [s for s in range(n) if s not in J]
            inp = ''.join([L[s] for s in range(n)] + [L[s].upper() if s not in J else L[s] for s in range(n)])
            out = ''.join([L[s] for s in keep] + [L[s].upper() for s in keep])
            R = np.einsum(inp + '->' + out, P)
            tot += ((x - 1) / 3) ** k * (np.abs(R) ** 2).sum()
    return tot / 2
def G_formula(n, x):
    B = 2 * 6 ** n / (3 ** n + 2 ** n) ** 2; R = (x + 2) / 3; S = (x - 1) / 3
    return 0.5 * ((2 - B) * R ** n + B + (2 + B) * S ** n)
err = max(abs(G_direct(n, x) - G_formula(n, x)) for n in range(2, 7) for x in (-0.5, 0.3, 0.2 + 0.7j, -1.3j))
print(f"(A) formula (89.1) versus partial traces, n = 2..6: max difference {err:.1e}")
def newton(f, x, h=1e-7, it=60):
    for _ in range(it):
        d = (f(x + h) - f(x - h)) / (2 * h); x = x - f(x) / d
    return x
n = 20
x20 = newton(lambda x: G_formula(n, x), -1 / 40 + 0.313j)
print(f"(A) n = 20: zero x = {x20:.6f}, |x - (-1/40 + 0.313i)| = {abs(x20 - (-1/40 + 0.313j)):.2e}, |x| = {abs(x20):.4f}")
for n in (40, 80, 160):
    xn = newton(lambda x: G_formula(n, x), 2j * np.pi / n)
    print(f"(A) n = {n}: zero x_n with n x_n = {n * xn:.5f} (2 pi i = {2j * np.pi:.5f})")
ok = all((Fr(2) ** (1 - n) + Fr(6 ** n, (3 ** n + 2 ** n) ** 2) if n % 2 == 0 else
          Fr(2 * 6 ** n, (3 ** n + 2 ** n) ** 2) * (Fr(1, 2) - Fr(1, 2 ** n))) > 0 for n in range(2, 200))
print(f"(A) endpoint values G(-1/2) > 0 for n = 2..199 (exact): {ok};"
      f" formula vs direct at n = 5: {G_formula(5, -0.5):.12f} vs {G_direct(5, -0.5):.12f}")
# ---------------- (B)
z = np.zeros((3, 3, 3))
for (a, b, c), s in {(0, 0, 0): 1, (1, 0, 1): 1, (2, 0, 2): 1, (0, 1, 1): -1, (1, 1, 0): 1, (0, 2, 2): -1, (2, 2, 0): 1}.items():
    z[a, b, c] = s
rho = np.einsum('abc,dbe->acde', z, z).reshape(9, 9)                 # Tr_B zz^T on A (x) C
rhoTA = rho.reshape(3, 3, 3, 3).transpose(2, 1, 0, 3).reshape(9, 9)
T = np.trace(rho @ rhoTA)
pur = [np.trace(np.einsum(sub, z, z) @ np.einsum(sub, z, z)) for sub in ('abc,dbc->ad', 'abc,adc->bd', 'abc,abd->cd')]
spec = [np.round(np.linalg.eigvalsh(np.einsum(sub, z, z)), 6).tolist() for sub in ('abc,dbc->ad', 'abc,adc->bd', 'abc,abd->cd')]
print(f"(B) ||z||^2 = {(z ** 2).sum():.0f}, T = Tr(rho rho^TA) = {T:.6f}, marginal purities {np.round(pur, 6).tolist()}, spectra {spec}")
def L3(X):
    return 3 * np.trace(X) * np.eye(3) - X - X.T
def apply_site(f, X, nn, s):
    Tt = X.reshape([3] * (2 * nn)); perm = [s, nn + s] + [i for i in range(2 * nn) if i not in (s, nn + s)]
    Tp = np.transpose(Tt, perm).reshape(3, 3, -1)
    out = np.stack([f(Tp[:, :, kk]) for kk in range(Tp.shape[2])], axis=2)
    return np.transpose(out.reshape([3, 3] + [3] * (2 * nn - 2)), np.argsort(perm)).reshape(3 ** nn, 3 ** nn)
zv = z.reshape(-1); X = np.outer(zv, zv)
for s in range(3): X = apply_site(L3, X, 3, s)
print(f"(B) critical scalar z^T L3^(x3)(zz^T) z = {zv @ X @ zv:.6f} (audit: 493)")
def overlap_ratio(R, d, k):
    RT = R.reshape(d, k, d, k).transpose(2, 1, 0, 3).reshape(d * k, d * k)
    return np.trace(R @ RT) / np.trace(R @ R)
worst = min(overlap_ratio((lambda M: M @ M.T)(rng.normal(size=(12, int(rng.integers(1, 13))))), 3, 4) for _ in range(20000))
print(f"(B) bound (89.2), d = 3, k = 4: smallest ratio over 20000 random real PSD = {worst:.5f} >= 2/3 - 1 = {-1/3:.5f}")
qi = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]], float)
qj = np.array([[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0]], float)
U = [np.eye(4), qi, qj]
assert np.allclose(qi @ qj, -qj @ qi) and np.allclose(qi @ qi.T, np.eye(4))
Rq = np.block([[U[a] @ U[b].T for b in range(3)] for a in range(3)])
print(f"(B) quaternion equality case: purity {np.trace(Rq @ Rq):.1f}, overlap {np.trace(Rq @ Rq.reshape(3, 4, 3, 4).transpose(2, 1, 0, 3).reshape(12, 12)):.1f},"
      f" ratio {overlap_ratio(Rq, 3, 4):+.5f}")
# ---------------- (C)
def Bloc(i, k, t):
    e = np.eye(3); return np.outer(e[i], e[k]) + np.outer(e[k], e[i]) + t * e[i][k] * np.eye(3)
def Q(z, nn, abc):
    a, b, c = abc
    tau = np.zeros((3, 3, 3, 3))
    for k1, k2, k3, k4 in itertools.product(range(3), repeat=4):
        tau[k1, k2, k3, k4] = a * (k1 == k2) * (k3 == k4) + b * (k1 == k3) * (k2 == k4) + c * (k1 == k4) * (k2 == k3)
    Z = z.reshape([3] * nn); out = np.einsum(Z, list(range(nn)), Z, list(range(nn, 2 * nn)),
                                             Z, list(range(2 * nn, 3 * nn)), Z, list(range(3 * nn, 4 * nn)),
                                             *sum(([tau, [s, nn + s, 2 * nn + s, 3 * nn + s]] for s in range(nn)), []), [])
    return float(out)
maxerr = 0; ddvv_ok = True
for nn in (1, 2):
    for t in (-0.7, 0.0, 0.4, 1.3):
        cc = 3 * t ** 2 + 4 * t + 8
        h = (3 * t ** 4 + 8 * t ** 3 + 16 * t ** 2 + 8 * t + 13, 8 * t ** 2 + 20 * t, 8 * t ** 2 + 20 * t + 25)
        j = (8 * t ** 2 + 20 * t + 9, 3 * t ** 4 + 8 * t ** 3 + 16 * t ** 2 + 8 * t + 2, 8 * t ** 2 + 20 * t + 9)
        for _ in range(3):
            zz = rng.normal(size=3 ** nn)
            Bs = []
            for I in itertools.product(range(3), repeat=nn):
                M = 0
                for K in itertools.product(range(3), repeat=nn):
                    term = np.array([[1.0]])
                    for s in range(nn): term = np.kron(term, Bloc(I[s], K[s], t))
                    M = M + zz[np.ravel_multi_index(K, [3] * nn)] * term
                Bs.append(M)
            S1 = sum((Bm ** 2).sum() for Bm in Bs); Msq = sum(Bm @ Bm for Bm in Bs)
            S2 = np.trace(Msq @ Msq); S3 = sum(np.trace(Bi @ Bj @ Bi @ Bj) for Bi in Bs for Bj in Bs)
            maxerr = max(maxerr, abs(S1 - cc ** nn * (zz ** 2).sum()), abs(S2 - Q(zz, nn, h)), abs(S3 - Q(zz, nn, j)))
            comm = sum(((Bi @ Bj - Bj @ Bi) ** 2).sum() for Bi in Bs for Bj in Bs)
            ddvv_ok &= comm <= S1 ** 2 + 1e-8
print(f"(C) DDVV identities (n = 1, 2; four t; random z): max error {maxerr:.1e}; DDVV inequality holds: {ddvv_ok}")
roots = np.roots([3, 8, 40, 68, 29])
for t in roots[np.abs(roots.imag) < 1e-9].real:
    j = np.array([8 * t ** 2 + 20 * t + 9, 3 * t ** 4 + 8 * t ** 3 + 16 * t ** 2 + 8 * t + 2, 8 * t ** 2 + 20 * t + 9])
    cc = 3 * t ** 2 + 4 * t + 8
    print(f"(C) root t = {t:.6f}: j = {np.round(j, 6).tolist()} -> j/(-j_1) = {np.round(j / -j[0], 6).tolist()},"
          f" error base c^2/kappa = {cc ** 2 / -j[0]:.4f}")
