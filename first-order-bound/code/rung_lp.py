"""Linear-programming test of the graded levels (see the Remarks): can a level A_k >= C(n-1,k) 2^-k S_0 be written as a
nonnegative combination of known inequalities?  Inputs, summed over all placements: the two- and three-copy theorems
for groupings of the copies (other copies spectators), the bounds S_m <= [C(n-1,m)+2C(n-1,m-1)] S_0, S_n <= 2 S_0, and,
with an extra variable N1 = ||C||_1^2 (S_0 <= N1 <= 2 S_0) and the spectator-free inequalities that
use it: the sharp two-factor bound for every full cut (K | K^c): ||Tr_K C||^2 + ||Tr_{K^c} C||^2 <= S_0 + N1/2 + S_n/2,
and Theorem 1 of the paper: S_1 <= (n-1) S_0 + N1.  Variables ordered (S_0..S_n, N1)."""
import sys, itertools
from math import comb, factorial
import numpy as np
from scipy.optimize import linprog
n = int(sys.argv[1])
def group_ineq(sizes):
    r = len(sizes); tot = sum(sizes)
    places = factorial(n) // (factorial(n - tot) * int(np.prod([factorial(s) for s in sizes])))
    v = np.zeros(n + 1)
    for J in itertools.product([0, 1], repeat=r):
        m = sum(s for s, j in zip(sizes, J) if j)
        v[m] += (-0.5) ** sum(J) * places / comb(n, m)
    return v
ineqs, names = [], []
def pad(v): return np.append(v, 0.0)
for k in range(1, n):
    for l in range(k, n - k + 1):
        ineqs.append(pad(group_ineq([k, l]))); names.append(f"G2({k},{l})")
for k1 in range(1, n + 1):
    for k2 in range(k1, n + 1):
        for k3 in range(k2, n + 1):
            if k1 + k2 + k3 <= n:
                ineqs.append(pad(group_ineq([k1, k2, k3]))); names.append(f"G3({k1},{k2},{k3})")
for k in range(1, n):   # full cuts with flatness term, summed over K of size k
    v = np.zeros(n + 2); v[0] = comb(n, k); v[n + 1] = comb(n, k) / 2; v[n] += comb(n, k) / 2
    v[k] -= 1; v[n - k] -= 1      # sum over K of ||Tr_K C||^2 = S_k, and of ||Tr_{K^c} C||^2 = S_{n-k}
    ineqs.append(v); names.append(f"CUT+({k})")
sp = np.zeros(n + 2); sp[0] = n - 1; sp[n + 1] = 1; sp[1] = -1; ineqs.append(sp); names.append("spade")
v = np.zeros(n + 2); v[0] = 2; v[n + 1] = -1; ineqs.append(v); names.append("N1<=2S0")
v = np.zeros(n + 2); v[n + 1] = 1; v[0] = -1; ineqs.append(v); names.append("N1>=S0")
for m in range(2, n + 1):
    v = np.zeros(n + 2); v[0] = comb(n - 1, m) + 2 * comb(n - 1, m - 1); v[m] = -1; ineqs.append(v); names.append(f"CS{m}")
v = np.zeros(n + 2); v[0] = 2; v[n] = -1; ineqs.append(v); names.append("Sn<=2S0")
for m in range(n + 1):
    v = np.zeros(n + 2); v[m] = 1; ineqs.append(v); names.append(f"S{m}>=0")
G = np.array(ineqs).T
for k in range(2, n + 1):
    t = np.zeros(n + 2)
    for m in range(k + 1): t[m] = (-0.5) ** m * comb(n - m, k - m)
    t[0] -= comb(n - 1, k) * 2.0 ** (-k)
    res = linprog(np.zeros(G.shape[1]), A_eq=G, b_eq=t, bounds=[(0, None)] * G.shape[1], method='highs')
    used = [(names[j], round(res.x[j], 4)) for j in range(len(names)) if res.status == 0 and res.x[j] > 1e-9 and not names[j].startswith('S')]
    print(f"n={n} rung k={k}{' (endpoint)' if k == n else ''}: {'DERIVABLE ' + str(used) if res.status == 0 else 'not derivable'}")
