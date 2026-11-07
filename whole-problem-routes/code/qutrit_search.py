"""Search used for n = 6 in Observation 7.4 (the mode `search` of map_thresholds.py does n = 2..5 the same way).
Does the qutrit map P_a(X) = a Tr(X) I - (X + X^T)/4 have the same n-copy positivity threshold as its
real-qubit compression, a_n = (1 + 2^(1 - 1/n))/4?  For a in {a_n - 0.004, a_n - 0.001, a_n + 0.001, a_n + 0.004} the
least value of <y|P_a^n(xx^T)|y> over real unit x, y in (R^3)^n is searched by alternating bottom eigenvectors from:
the embedded GHZ_Y vector (|+i>^n + |-i>^n)/sqrt2 (real), the same after random rotations of each copy in R^3 and random
perturbations of size 0.05 and 0.3 in all three dimensions, and plain random starts; then continuation: the best pair
at each a seeds the next a.  Prints the least value found and the one-copy ranks of the best x.
usage: python qutrit_search.py n random_starts   (logs/qutrit_search_n6.log: python qutrit_search.py 6 4)"""
import sys
import numpy as np

n = int(sys.argv[1]); rs = int(sys.argv[2])
rng = np.random.default_rng(31337 + n)
d = 3; D = d ** n


def phi_for(a):
    e = np.eye(d); M = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            E = np.outer(e[i], e[j])
            M[:, d * i + j] = (a * np.trace(E) * np.eye(d) - (E + E.T) / 4).reshape(d * d)
    M4 = M.reshape(d, d, d, d)

    def Phi(X):
        T = X.reshape([d] * (2 * n))
        for c in range(n):
            T = np.moveaxis(np.tensordot(M4, T, axes=([2, 3], [c, n + c])), [0, 1], [c, n + c])
        return T.reshape(D, D)
    return Phi


def ghzy():
    p = np.array([1, 1j, 0]) / np.sqrt(2)
    g = np.ones(1, dtype=complex)
    for _ in range(n):
        g = np.kron(g, p)
    v = np.real(g + np.conj(g)); return v / np.linalg.norm(v)


def rot_all(v):
    T = v.reshape([d] * n)
    for c in range(n):
        Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
        T = np.moveaxis(np.tensordot(Q, T, axes=([1], [c])), 0, c)
    return T.reshape(D)


def descend(Phi, x):
    x = x / np.linalg.norm(x); val = np.inf; y = None
    for it in range(800):
        try:
            w, V = np.linalg.eigh(Phi(np.outer(x, x))); y = V[:, 0]
            w2, V2 = np.linalg.eigh(Phi(np.outer(y, y)))
        except np.linalg.LinAlgError:
            break
        nv = w2[0]; x = V2[:, 0]
        if abs(nv - val) < 1e-15:
            val = nv; break
        val = nv
    return val, x, y


def ranks(v):
    T = v.reshape([d] * n); out = []
    for c in range(n):
        A = np.moveaxis(T, c, 0).reshape(d, -1); ev = np.linalg.eigvalsh(A @ A.T)
        out.append(int(np.sum(ev > 1e-7)))
    return out


an = (1 + 2 ** (1 - 1 / n)) / 4
g = ghzy()
carry = None
for a in (an - 0.004, an - 0.001, an + 0.001, an + 0.004):
    Phi = phi_for(a)
    seeds = [g, rot_all(g)] + [rot_all(g) + eps * rng.normal(size=D) for eps in (0.05, 0.3) for _ in range(4)]
    seeds += [rng.normal(size=D) for _ in range(rs)]
    if carry is not None:
        seeds = [carry, carry + 0.05 * rng.normal(size=D)] + seeds
    best = (np.inf, None)
    for s in seeds:
        v, x, y = descend(Phi, s)
        if v < best[0]:
            best = (v, x)
    carry = best[1]
    m = a - 0.5

    print(f'n = {n}, a = {a:.6f} (qubit threshold {an:.6f}): least value found {best[0]:+.4e} (m^n = {m ** n:.3e}); '
          f'one-copy ranks of x {ranks(best[1])}', flush=True)
