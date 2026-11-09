"""Exact checks for the identity-count polytope (Section 7, "The identity-count polytope").
pibar_C(k) = sum_{|K|=k} ||C_K||^2/||C||^2, C_K the part of C that is the identity exactly on the copies in K.
(1) Realise the four vertices of S_2 and the identified vertices of S_3 by explicit operators, in exact arithmetic
    (Fractions), from the formula pibar via Moebius inversion of the partial-trace norms.
(2) Check the generating-function factorisations G(x) = sum_k pibar(k) x^k of the zero vertices.
usage: python check_vertices.py"""
import itertools
from fractions import Fraction as Fr
import numpy as np
def pibar_exact(C, N):
    """C: dict {(row, col): Fraction or complex-rational as (re, im) Fractions} -> exact pibar; uses real C only."""
    D = 3 ** N
    M = np.array([[C.get((i, j), Fr(0)) for j in range(D)] for i in range(D)], dtype=object)
    subsets = list(itertools.product([0, 1], repeat=N))
    def tr_norm(K):
        idx = list(itertools.product(range(3), repeat=N))
        kept = [s for s in range(N) if not K[s]]
        acc = {}
        for a, ia in zip(idx, range(D)):
            for b, ib in zip(idx, range(D)):
                if M[ia, ib] == 0: continue
                if any(a[s] != b[s] for s in range(N) if K[s]): continue
                key = (tuple(a[s] for s in kept), tuple(b[s] for s in kept))
                acc[key] = acc.get(key, Fr(0)) + M[ia, ib]
        return sum(v * v for v in acc.values())
    mu = {K: tr_norm(K) / Fr(3) ** sum(K) for K in subsets}
    e = {L: sum(((-1) ** (sum(K) - sum(L))) * mu[K] for K in subsets if all(K[s] >= L[s] for s in range(N))) for L in subsets}
    nrm = sum(M[i, j] * M[i, j] for i in range(D) for j in range(D))
    return [sum(e[L] for L in subsets if sum(L) == k) / nrm for k in range(N + 1)]
def idx(*d): return int(''.join(map(str, d)), 3)
def show(name, C, N):
    p = pibar_exact(C, N); den = 3 ** N
    print(f"  {name:48s} pibar x{den} = {[str(x * den) for x in p]}   endpoint = {sum(x * Fr(-1, 2) ** k for k, x in enumerate(p))}")
print("Two copies (vertices of S_2):")
show("|0><1| (x) |0><1|  (traceless)", {(idx(0, 0), idx(1, 1)): Fr(1)}, 2)
show("Pi_2 (x) |0><1|", {(idx(0, 0), idx(0, 1)): Fr(1), (idx(1, 0), idx(1, 1)): Fr(1)}, 2)
show("Pi_2 (x) |0><0|", {(idx(0, 0), idx(0, 0)): Fr(1), (idx(1, 0), idx(1, 0)): Fr(1)}, 2)
# balanced rank-two projector: span{(|00>+|11>+|22>)/sqrt3, (|01>+|12>+|20>)/sqrt3}  (real, entries 1/3)
B = {}
for v in [[(0, 0), (1, 1), (2, 2)], [(0, 1), (1, 2), (2, 0)]]:
    for a in v:
        for b in v: B[(idx(*a), idx(*b))] = B.get((idx(*a), idx(*b)), Fr(0)) + Fr(1, 3)
show("balanced projector span{Phi, X-shifted Phi}", B, 2)
print("Three copies:")
show("Pi_2 (x) |Phi><Phi|  (zero vertex)", {(idx(c, a, a), idx(c, b, b)): Fr(1, 3) for c in range(2) for a in range(3) for b in range(3)}, 3)
show("|000><000| + |111><111|  (GHZ zero)", {(idx(0, 0, 0), idx(0, 0, 0)): Fr(1), (idx(1, 1, 1), idx(1, 1, 1)): Fr(1)}, 3)
show("|000><000| - |111><111|", {(idx(0, 0, 0), idx(0, 0, 0)): Fr(1), (idx(1, 1, 1), idx(1, 1, 1)): Fr(-1)}, 3)
show("Pi_2 (x) |0><1| (x) |0><1|", {(idx(c, 0, 0), idx(c, 1, 1)): Fr(1) for c in range(2)}, 3)
show("Pi_2 (x) |0><0| (x) |0><1|", {(idx(c, 0, 0), idx(c, 0, 1)): Fr(1) for c in range(2)}, 3)
show("Pi_2 (x) |0><0| (x) |0><0|", {(idx(c, 0, 0), idx(c, 0, 0)): Fr(1) for c in range(2)}, 3)
B3 = {}
for a in B:
    for c in [(0, 1)]:
        B3[(a[0] * 3 + 0, a[1] * 3 + 1)] = B[a]
show("balanced 2-copy projector (x) |0><1|", B3, 3)
