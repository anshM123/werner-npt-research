"""Section 5: activation of two-pair Werner mixtures by one more Werner copy.

   sigma_c = sum_S c_S Pi_S on two qutrit pairs, Pi_S = products of the (unnormalised) symmetric / antisymmetric
   projectors; rho_W(1/2)^{(x)2} is c = (c_ss, c_sa, c_as, c_aa) = (1, 3, 3, 9). The form of sigma_c (x) rho_W is
   Q(C) = sum_J gamma_J ||Tr_J C||^2 with gamma_{J1 J2 J3} = gamma^c_{J1 J2} (-1/2)^{J3}, where
   gamma^c_J = 2^{-2} sum_S c_S (-1)^{|J cap S|}. The script checks the closed form
       Q(C0)/||C0||^2 = c_ss/2 - c_sa/6 + 3 c_as/4 - c_aa/4
   for C0 = |Phi><Phi| on (pair 1, extra copy) (x) P on pair 2, and evaluates the two-pair minimum (1-copy
   2-undistillability of sigma_c) by local searches at a few points of the slice c = (1, x, x, y)."""
import itertools
import numpy as np
from scipy.optimize import minimize

d = 3
rng = np.random.default_rng(6)


def gam2(c):
    return {J: sum(cS * (-1) ** sum(a * b for a, b in zip(J, S)) for S, cS in c.items()) / 4
            for J in itertools.product([0, 1], repeat=2)}


def form(X, k, coef):
    T = X.reshape([d] * (2 * k)); val = 0.0
    for J, cJ in coef.items():
        A, mm = T, k
        for i in reversed(range(k)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + mm); mm -= 1
        val += cJ * np.vdot(A, A).real
    return val


def with_werner(g):
    return {J + (j3,): gJ * (-0.5) ** j3 for J, gJ in g.items() for j3 in (0, 1)}


Phi = np.eye(d).reshape(-1) / np.sqrt(d)                       # on (pair 1, extra copy)
P = np.diag([1.0, 1.0, 0.0])
C0 = np.einsum('ac,bd,ef->aecbfd', Phi.reshape(d, d), Phi.reshape(d, d).conj(), P).reshape(27, 27)
# index check: C0[(i0 i1 i2),(j0 j1 j2)] = Phi(i0, i2) conj Phi(j0, j2) P[i1, j1]
for c in [(1, 3, 3, 9), (1, 2, 2, 7), (1, 1, 1, 5), (1, 0.7, 2.3, 4.1), (0.4, 1.9, 0.6, 3.3)]:
    cd = dict(zip([(0, 0), (0, 1), (1, 0), (1, 1)], c))
    direct = form(C0, 3, with_werner(gam2(cd))) / np.vdot(C0, C0).real
    closed = c[0] / 2 - c[1] / 6 + 3 * c[2] / 4 - c[3] / 4
    print(f"c = {c}: direct {direct:+.10f}   closed form {closed:+.10f}", flush=True)


def min_rank2(coef, k, starts=16):
    D = d ** k; best = np.inf
    for _ in range(starts):
        def f(x):
            Y = (x[:2 * D] + 1j * x[2 * D:4 * D]).reshape(D, 2); Z = (x[4 * D:6 * D] + 1j * x[6 * D:]).reshape(D, 2)
            C = Y @ Z.conj().T
            return form(C, k, coef) / np.vdot(C, C).real
        r = minimize(f, rng.standard_normal(8 * D), method='L-BFGS-B', options={'maxiter': 4000})
        best = min(best, r.fun)
    return best


print("\nslice c = (1, x, x, y): two-pair minimum m2 (local searches) and the witness value 1/2 + 7x/12 - y/4")
for x, y in [(2, 6.75), (2, 20 / 3), (1, 4.5), (3, 9)]:
    cd = {(0, 0): 1, (0, 1): x, (1, 0): x, (1, 1): y}
    print(f"  x={x}, y={y:.4f}: m2 = {min_rank2(gam2(cd), 2):+.5f}   witness {0.5 + 7 * x / 12 - y / 4:+.5f}", flush=True)
