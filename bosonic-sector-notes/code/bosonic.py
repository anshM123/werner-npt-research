"""Shared helpers for the bosonic sector, End(Sym^n C^3), in the normalised occupation basis.

    q_n(t, C) = sum_J (-t)^{|J|} ||Tr_J C||_2^2
              = sum_j (-t)^j (n-j)!/(j! n!) ||A^j C||^2,    A(C) = sum_i a_i C a_i^dag,

where a_i : Sym^m -> Sym^{m-1} are the annihilation operators."""
import itertools
import numpy as np
from math import factorial, comb
from fractions import Fraction


def basis(n):
    """occupation vectors (n0, n1, n2) with n0 + n1 + n2 = n"""
    return [(a, b, n - a - b) for a in range(n + 1) for b in range(n + 1 - a)]


_ANN = {}


def ann(m):
    """annihilation matrices a_0, a_1, a_2 : Sym^m -> Sym^{m-1}"""
    if m not in _ANN:
        B, Bm = basis(m), basis(m - 1)
        idx = {t: k for k, t in enumerate(Bm)}
        out = []
        for i in range(3):
            A = np.zeros((len(Bm), len(B)))
            for k, t in enumerate(B):
                if t[i] > 0:
                    tt = list(t); tt[i] -= 1
                    A[idx[tuple(tt)], k] = np.sqrt(t[i])
            out.append(A)
        _ANN[m] = out
    return _ANN[m]


def A_op(X, m):
    """A : End(Sym^m) -> End(Sym^{m-1})"""
    return sum(a @ X @ a.T for a in ann(m))


def A_adj(X, m):
    """adjoint of A : End(Sym^{m-1}) -> End(Sym^m)"""
    return sum(a.T @ X @ a for a in ann(m))


def coeff(n, j, t):
    return (-t) ** j * factorial(n - j) / (factorial(j) * factorial(n))


def q_and_grad(C, n, t=0.5):
    """q_n(t, C) and its Wirtinger gradient Q(C) = dq/d(conj C)"""
    Xs = [C]
    for j in range(1, n + 1):
        Xs.append(A_op(Xs[-1], n - j + 1))
    q = sum(coeff(n, j, t) * np.vdot(X, X).real for j, X in enumerate(Xs))
    G = coeff(n, n, t) * Xs[n]
    for j in range(n, 0, -1):
        G = coeff(n, j - 1, t) * Xs[j - 1] + A_adj(G, n - j + 1)
    return q, G


def gamma(n, k, t=Fraction(1, 2)):
    """eigenvalue of the q-form on the SU(3) component V_(k,k) of End(Sym^n C^3)"""
    return sum((-t) ** j * comb(n - k, j) * comb(n + k + 2, j) / Fraction(comb(n, j)) for j in range(n - k + 1))


def sym_embed(n):
    """isometry Sym^n -> (C^3)^{(x)n}"""
    B = basis(n)
    E = np.zeros((3 ** n, len(B)))
    for k, m in enumerate(B):
        word = [0] * m[0] + [1] * m[1] + [2] * m[2]
        perms = set(itertools.permutations(word))
        for p in perms:
            E[int(''.join(map(str, p)), 3), k] = 1 / np.sqrt(len(perms))
    return E


def q_direct(X, n, t=0.5, d=3):
    """q_n(t, X) from the definition, X an operator on (C^d)^{(x)n}"""
    T = X.reshape([d] * (2 * n))
    tot = 0.0
    for J in itertools.product([0, 1], repeat=n):
        A, mm = T, n
        for i in reversed(range(n)):
            if J[i]:
                A = np.trace(A, axis1=i, axis2=i + mm); mm -= 1
        M = A.reshape(d ** mm, d ** mm)
        tot += (-t) ** sum(J) * np.vdot(M, M).real
    return tot


def min_rank2(n, t=0.5, starts=12, mask=None, rng=None):
    """minimum of q_n(t,C)/||C||^2 over C = Y Z^*, Y, Z : C^2 -> Sym^n (rows restricted by mask if given)"""
    from scipy.optimize import minimize
    rng = rng if rng is not None else np.random.default_rng(0)
    D = len(basis(n))
    m = np.ones(D) if mask is None else mask.astype(float)

    def f(x):
        Y = (x[:2 * D] + 1j * x[2 * D:4 * D]).reshape(D, 2) * m[:, None]
        Z = (x[4 * D:6 * D] + 1j * x[6 * D:]).reshape(D, 2) * m[:, None]
        C = Y @ Z.conj().T
        nn = np.vdot(C, C).real
        q, G = q_and_grad(C, n, t)
        val = q / nn
        Gc = (G - val * C) / nn
        GY = (Gc @ Z) * m[:, None]
        GZ = (Gc.conj().T @ Y) * m[:, None]
        return val, 2 * np.concatenate([GY.real.ravel(), GY.imag.ravel(), GZ.real.ravel(), GZ.imag.ravel()])

    return min(minimize(f, rng.standard_normal(8 * D), jac=True, method='L-BFGS-B',
                        options={'maxiter': 20000, 'gtol': 1e-12, 'ftol': 1e-15}).fun for _ in range(starts))
