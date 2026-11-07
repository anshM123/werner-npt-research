"""Exploration: search for small values of <y|P_a^{(x)n}(x x^T)|y> over real open-boundary MPS x, y (bond dimension D)
by alternating one-site DMRG: for fixed y, M_y = P_a^{(x)n}(y y^T) is an MPO (bond D^2) and x is updated site by site as the
lowest eigenvector of the effective 3D^2-dimensional problem (mixed canonical form, so the norm is local); then the roles
of x and y swap (P_a is self-adjoint, so <y|P(xx^T)|y> = <x|P(yy^T)|x>).
usage: python dmrg_search.py n D starts sweeps a1 a2 ...    (prints the least value found for each a, with the product
value (a-1/2)^n and the witness value 2^-n[(2a-1/2)^n - 1/2] of Theorem 7.1)"""
import sys, time
import numpy as np

n, D, starts, sweeps = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
avals = [float(v) for v in sys.argv[5:]]
rng = np.random.default_rng(7000 + 17 * n + D)


def Pten(a):     # P[k, l, i, j] = <k| P_a(e_i e_j^T) |l>
    P = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    P[k, l, i, j] = a * (i == j) * (k == l) - ((k == i) * (l == j) + (k == j) * (l == i)) / 4
    return P


def dims(n, D):
    d = [1]
    for s in range(1, n):
        d.append(min(D, 3 ** s, 3 ** (n - s)))
    d.append(1)
    return d


def random_mps(n, D, thick=False):
    d = dims(n, D)
    A = [rng.normal(size=(d[s], 3, d[s + 1])) for s in range(n)]
    if thick:
        for s in range(n):
            A[s] += 0.5 * rng.normal(size=A[s].shape)
    return A


def ghz_mps(n, D):
    # real GHZ_Y-type state (p^n + conj(p)^n)/sqrt2 written as a bond-2 real MPS, padded
    d = dims(n, D)
    c, s_ = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    A = []
    for k in range(n):
        T = np.zeros((d[k], 3, d[k + 1]))
        # rotation-angle representation: p^{(x)n} ~ (e1 + i e2)^n, real part as 2x2 transfer
        R = np.zeros((2, 3, 2)); R[0, :, 0] = c; R[1, :, 1] = c; R[0, :, 1] = -s_; R[1, :, 0] = s_
        if k == 0:
            T[0, :, :min(2, d[1])] = R[0, :, :min(2, d[1])]
        elif k == n - 1:
            T[:min(2, d[k]), :, 0] = R[:min(2, d[k]), :, 0]
        else:
            T[:min(2, d[k]), :, :min(2, d[k + 1])] = R[:min(2, d[k]), :, :min(2, d[k + 1])]
        A.append(T + 1e-3 * rng.normal(size=T.shape))
    return A


def mpo_from(B, P):
    # W_s[a, b, i, j] with a, b MPO bond indices (pairs of B bonds): sum_{k,l} P[k,l,i,j] B[k] (x) B[l]
    W = []
    for T in B:
        dl, _, dr = T.shape
        Ws = np.einsum('klij,akb,cld->acbdij', P, T, T).reshape(dl * dl, dr * dr, 3, 3)
        W.append(Ws)
    return W


def left_canonical_sweep(A):
    A = [t.copy() for t in A]
    for s in range(len(A) - 1):
        dl, p, dr = A[s].shape
        Q, R = np.linalg.qr(A[s].reshape(dl * p, dr))
        k = Q.shape[1]
        A[s] = Q.reshape(dl, p, k)
        A[s + 1] = np.einsum('ab,bpc->apc', R, A[s + 1])
    A[-1] /= np.linalg.norm(A[-1])
    return A


def dmrg_min(W, A, nsweep):
    """lowest value of <x|M|x>/<x|x> over MPS x with the bond dimensions of A (one-site sweeps)."""
    n = len(A)
    A = left_canonical_sweep(A)
    # right-canonicalise by sweeping from the right
    for s in range(n - 1, 0, -1):
        dl, p, dr = A[s].shape
        Q, R = np.linalg.qr(A[s].reshape(dl, p * dr).T)
        k = Q.shape[1]
        A[s] = Q.T.reshape(k, p, dr)
        A[s - 1] = np.einsum('apb,cb->apc', A[s - 1], R)
    A[0] /= np.linalg.norm(A[0])
    Lenv = [None] * (n + 1); Renv = [None] * (n + 1)
    Lenv[0] = np.ones((1, 1, 1)); Renv[n] = np.ones((1, 1, 1))
    for s in range(n - 1, -1, -1):
        Renv[s] = np.einsum('apb,cdpq,eqf,bdf->ace', A[s], W[s], A[s], Renv[s + 1])
    val = None
    for sw in range(nsweep):
        for direction in (1, -1):
            rng_sites = range(n) if direction == 1 else range(n - 1, -1, -1)
            for s in rng_sites:
                dl, p, dr = A[s].shape
                H = np.einsum('ace,cdpq,bdf->apbeqf', Lenv[s], W[s], Renv[s + 1]).reshape(dl * p * dr, dl * p * dr)
                H = (H + H.T) / 2
                w, V = np.linalg.eigh(H)
                val = w[0]
                T = V[:, 0].reshape(dl, p, dr)
                if direction == 1 and s < n - 1:
                    Q, R = np.linalg.qr(T.reshape(dl * p, dr))
                    A[s] = Q.reshape(dl, p, Q.shape[1])
                    A[s + 1] = np.einsum('ab,bpc->apc', R, A[s + 1])
                    Lenv[s + 1] = np.einsum('apb,cdpq,eqf,ace->bdf', A[s], W[s], A[s], Lenv[s])
                elif direction == -1 and s > 0:
                    Q, R = np.linalg.qr(T.reshape(dl, p * dr).T)
                    A[s] = Q.T.reshape(Q.shape[1], p, dr)
                    A[s - 1] = np.einsum('apb,cb->apc', A[s - 1], R)
                    Renv[s] = np.einsum('apb,cdpq,eqf,bdf->ace', A[s], W[s], A[s], Renv[s + 1])
                else:
                    A[s] = T
    return val, A


T0 = time.time()
for a in avals:
    P = Pten(a)
    best = np.inf
    seeds = [('ghz', ghz_mps(n, D))] + [('rand', random_mps(n, D)) for _ in range(starts)] + \
            [('thick', random_mps(n, D, thick=True)) for _ in range(max(1, starts // 2))]
    for kind, x in seeds:
        y = [t.copy() for t in x] if kind == 'ghz' else random_mps(n, D)
        val = np.inf
        for it in range(sweeps):
            v1, y = dmrg_min(mpo_from(x, P), y, 2)      # best y for this x
            v2, x = dmrg_min(mpo_from(y, P), x, 2)      # best x for this y
            if abs(v2 - val) < 1e-13 * max(1, abs(val)):
                val = v2; break
            val = v2
        best = min(best, val)
    an = (1 + 2 ** (1 - 1 / n)) / 4
    print(f'n = {n}, D = {D}, a = {a:.5f} (a_n = {an:.5f}): least value {best:+.4e}; product {(a - 0.5) ** n:.4e}; '
          f'witness {2.0 ** -n * ((2 * a - 0.5) ** n - 0.5):+.4e}  [{time.time() - T0:.0f}s]', flush=True)
