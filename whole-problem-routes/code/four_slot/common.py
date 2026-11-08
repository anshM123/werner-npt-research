"""Shared helpers: four-slot (two-replica-squared) operators for the crossed Cauchy-Schwarz form.

Slots per site: 0 = a, 1 = d, 2 = c, 3 = b   (Phi = a (x) d (x) c (x) b).
G = W_t^{(x)n},  W_t = I - t F  on (replica1, replica2).
  <Phi| A |Phi> = <ab|G|ab><cd|G|cd>,     A_1 = W_{03} W_{12}
  <Phi| B |Phi> = |<ad|G|cb>|^2,           B_1 = W_{01} W_{23} P_{(02)(13)}
Crossed CS at n copies  <=>  A_1^{(x)n} - B_1^{(x)n} >= 0 on 4-party product vectors.
"""
import itertools
import numpy as np

d = 3
NS = 4  # slots


def perm_op(perm, d=3, ns=4):
    """Operator P with P |i_0..i_{ns-1}> = |j> where j_{perm[k]} = i_k (slot k moved to slot perm[k])."""
    D = d ** ns
    P = np.zeros((D, D))
    for idx in itertools.product(range(d), repeat=ns):
        j = [0] * ns
        for k in range(ns):
            j[perm[k]] = idx[k]
        a = np.ravel_multi_index(idx, (d,) * ns)
        b = np.ravel_multi_index(tuple(j), (d,) * ns)
        P[b, a] = 1.0
    return P


def swap(i, j, d=3, ns=4):
    p = list(range(ns))
    p[i], p[j] = j, i
    return perm_op(p, d, ns)


def ptrans(X, S, d=3, ns=4):
    """Partial transpose on the slots in S."""
    if len(S) == 0:
        return X.copy()
    T = X.reshape((d,) * (2 * ns))
    axes = list(range(2 * ns))
    for k in S:
        axes[k], axes[ns + k] = axes[ns + k], axes[k]
    return T.transpose(axes).reshape(d ** ns, d ** ns)


def AB(t, d=3):
    I = np.eye(d ** 4)
    W = lambda i, j: I - t * swap(i, j, d)
    A = W(0, 3) @ W(1, 2)
    B = W(0, 1) @ W(2, 3) @ perm_op([2, 3, 0, 1], d)
    return A, B


def all_perms(d=3, ns=4):
    perms = list(itertools.permutations(range(ns)))
    return perms, [perm_op(list(p), d, ns) for p in perms]


def inv_perm(p):
    q = [0] * len(p)
    for i, v in enumerate(p):
        q[v] = i
    return tuple(q)


def subsets(ns=4):
    out = []
    for r in range(ns + 1):
        out += [tuple(c) for c in itertools.combinations(range(ns), r)]
    return out


def rand_vec(dim, rng):
    v = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    return v / np.linalg.norm(v)
