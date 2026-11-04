"""Shared helpers: the endpoint form q_n by direct partial traces, and the index sets of the class."""
import itertools
import numpy as np

d = 3


def q_direct(X, n, t=0.5):
    """q_n(t, X) = sum_J (-t)^|J| ||Tr_J X||_F^2, computed as <X, Psi^{(x)n}(X)> with Psi(Y) = Y - t Tr(Y) I."""
    T = X.reshape([d] * (2 * n))
    for s in range(n):
        tr = np.trace(T, axis1=s, axis2=n + s)
        full = np.moveaxis(np.einsum('ab,...->...ab', np.eye(d), tr), [2 * n - 2, 2 * n - 1], [s, n + s])
        T = T - t * full
    return np.vdot(X, T.reshape(d ** n, d ** n)).real


def strings(n):
    return list(itertools.product(range(d), repeat=n))


def weight(st):
    return sum(v != 0 for v in st)


def single_excitations(n):
    """Rows of the class: (s, z) <-> the string with value z in {1,2} at site s and 0 elsewhere."""
    return [(s, z) for s in range(n) for z in (1, 2)]


def row_string(n, s, z):
    return tuple(z if r == s else 0 for r in range(n))


def sector_form(Xs, n, c, rows, cols):
    """Right side of the sector formula, for X with entries Xs[i, j] = X[rows[i]; cols[j]] (cols: strings of weight c):
    4 Off + 2 Diag + sum over rho of weight c-1 of (|T_rho|^2 - 2 sum_{s free} |tau_{rho,s}|^2)."""
    ridx = {r: i for i, r in enumerate(rows)}
    cidx = {st: j for j, st in enumerate(cols)}
    off = diag = 0.0
    for i, (s, z) in enumerate(rows):
        for j, sig in enumerate(cols):
            if sig[s] != 0:
                diag += abs(Xs[i, j]) ** 2
            else:
                off += abs(Xs[i, j]) ** 2
    tot = 4 * off + 2 * diag
    for rho in strings(n):
        if weight(rho) != c - 1:
            continue
        taus = []
        for s in range(n):
            if rho[s] != 0:
                continue
            tau = 0
            for z in (1, 2):
                sig = list(rho); sig[s] = z
                tau += Xs[ridx[(s, z)], cidx[tuple(sig)]]
            taus.append(tau)
        taus = np.array(taus)
        tot += abs(taus.sum()) ** 2 - 2 * np.sum(abs(taus) ** 2)
    return tot


def embed(Xs, n, rows, cols):
    """Place the class matrix Xs into a full operator on (C^3)^{(x)n}."""
    idx = {st: i for i, st in enumerate(strings(n))}
    X = np.zeros((d ** n, d ** n), complex)
    for i, (s, z) in enumerate(rows):
        for j, sig in enumerate(cols):
            X[idx[row_string(n, s, z)], idx[sig]] = Xs[i, j]
    return X
