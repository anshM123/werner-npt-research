"""Exact per-copy data for the three-copy thickness certificate.

For k = 2, 3, 4 strands of R^3: Brauer diagrams (perfect matchings of 2k points; top 0..k-1, bottom k..2k-1) as
integer matrices, exact integer bases V_j of the m = 0 vectors of spin j (kernel of the real generator L_z, eigenspaces
of the Casimir), the integer Gram-block map Phi[(j, r, c), a] = (V_j^T D_a V_j)[r, c], a right inverse Psi, and the
exact relations among diagrams (k = 4: 105 diagrams span 91 dimensions).  Results cached in exact_alg_k{k}.pkl."""
import itertools
import os
import pickle
from fractions import Fraction
import numpy as np
import math
import flint

D = 3


def matchings(points):
    if not points:
        yield []
        return
    p = points[0]
    for i in range(1, len(points)):
        rest = points[1:i] + points[i + 1:]
        for m in matchings(rest):
            yield [(p, points[i])] + m


def canon(m):
    return tuple(sorted(tuple(sorted(p)) for p in m))


def diagram_list(k):
    ident = canon([(s, s + k) for s in range(k)])
    return [ident] + [m for m in sorted(set(canon(m) for m in matchings(list(range(2 * k))))) if m != ident]


def diagram_matrix(m, k):
    M = np.zeros((D ** k, D ** k), dtype=np.int64)
    for ii in itertools.product(range(D), repeat=2 * k):
        if all(ii[p] == ii[q] for p, q in m):
            r = c = 0
            for s in range(k):
                r = r * D + ii[s]
                c = c * D + ii[k + s]
            M[r, c] = 1
    return M


def generators(k):
    eps = np.zeros((3, 3, 3), dtype=np.int64)
    for (i, j, l), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
        eps[i, j, l] = s
    L = [-eps[q] for q in range(3)]
    out = []
    for Lq in L:
        tot = np.zeros((D ** k, D ** k), dtype=np.int64)
        for s in range(k):
            f = [np.eye(3, dtype=np.int64)] * k
            f[s] = Lq
            mm = f[0]
            for x in f[1:]:
                mm = np.kron(mm, x)
            tot += mm
        out.append(tot)
    return out


def fz(M):
    return flint.fmpz_mat([[int(x) for x in row] for row in M])


def nullspace_int(M):
    """integer basis (as numpy object array, columns) of the null space of an integer matrix"""
    X, nul = fz(M).nullspace()
    cols = []
    for c in range(nul):
        col = [int(X[r, c]) for r in range(X.nrows())]
        g = 0
        for v in col:
            g = math.gcd(g, abs(v))
        col = [v // g for v in col]
        cols.append(col)
    return np.array(cols, dtype=object).T if cols else np.zeros((M.shape[1], 0), dtype=object)


def build(k):
    fn = f'exact_alg_k{k}.pkl'
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    diags = diagram_list(k)
    mats = [diagram_matrix(m, k) for m in diags]
    T = generators(k)
    N = nullspace_int(T[2])                                  # kernel of L_z (integer columns)
    Cas = T[0] @ T[0] + T[1] @ T[1] + T[2] @ T[2]
    Nq = flint.fmpq_mat(fz(N.astype(object)))
    NtN = Nq.transpose() * Nq
    kk = NtN.inv() * Nq.transpose() * flint.fmpq_mat(fz(Cas)) * Nq
    assert flint.fmpq_mat(fz(Cas)) * Nq == Nq * kk
    V = {}
    for j in range(k + 1):
        Mj = kk + j * (j + 1) * flint.fmpq_mat(kk.nrows(), kk.ncols(), [1 if r == c else 0 for r in range(kk.nrows()) for c in range(kk.ncols())])
        num, den = Mj.numer_denom()
        ns = nullspace_int(np.array([[int(num[r, c]) for c in range(num.ncols())] for r in range(num.nrows())], dtype=object))
        if ns.shape[1]:
            Vj = (np.array(N, dtype=object) @ ns)
            # make columns coprime integers
            for c in range(Vj.shape[1]):
                g = 0
                for v in Vj[:, c]:
                    g = math.gcd(g, abs(int(v)))
                Vj[:, c] = [int(v) // g for v in Vj[:, c]]
            V[j] = Vj.astype(np.int64)
    spins = sorted(V)
    mj = {j: V[j].shape[1] for j in spins}
    coords = [(j, r, c) for j in spins for r in range(mj[j]) for c in range(mj[j])]
    cidx = {c: i for i, c in enumerate(coords)}
    Phi = np.zeros((len(coords), len(diags)), dtype=object)
    for a, M in enumerate(mats):
        for j in spins:
            G = V[j].T @ M @ V[j]
            for r in range(mj[j]):
                for c in range(mj[j]):
                    Phi[cidx[(j, r, c)], a] = int(G[r, c])
    rank = fz(Phi).rank()
    assert rank == len(coords), (rank, len(coords))
    import scipy.linalg as sla
    _, _, piv = sla.qr(np.array(Phi, dtype=float), pivoting=True)
    basis = sorted(piv[:rank].tolist())
    Sub = flint.fmpq_mat(fz(Phi[:, basis]))
    Subinv = Sub.inv()
    Psi = [[Fraction(0)] * len(coords) for _ in range(len(diags))]
    for ii, a in enumerate(basis):
        for c in range(len(coords)):
            x = Subinv[ii, c]
            Psi[a][c] = Fraction(int(x.p), int(x.q))
    rels = nullspace_int(Phi)
    relations = [[int(x) for x in rels[:, c]] for c in range(rels.shape[1])]
    for r_ in relations:
        assert np.abs(sum(c * mats[a] for a, c in enumerate(r_))).max() == 0
    out = dict(k=k, diags=diags, mats=mats, V=V, spins=spins, mj=mj, coords=coords,
               Phi=Phi, Psi=Psi, basis=basis, relations=relations)
    pickle.dump(out, open(fn, 'wb'))
    return out


if __name__ == '__main__':
    import time
    for k in (2, 3, 4):
        t0 = time.time()
        A = build(k)
        dens = set()
        for row in A['Psi']:
            for x in row:
                dens.add(x.denominator)
        print(f'k = {k}: {len(A["diags"])} diagrams, spins {A["mj"]}, {len(A["coords"])} block coords, '
              f'{len(A["relations"])} relations, Psi denominators up to {max(dens)}  [{time.time() - t0:.1f}s]', flush=True)
