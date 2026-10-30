"""Exact feasibility of (rho^k)_k in the cone of symmetrised product generators at a rational rho
for products of per-site Werner operators X(r) = P^s + r P^a with at most m3 endpoint factors. Prints an exact certificate (basis + rational weights) if feasible.
usage: copy_window_check.py n d m3 rE_num rE_den rho_num rho_den
"""
import itertools
import sys
from fractions import Fraction as Fr
from math import comb, prod

n, d, m3 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rE = Fr(int(sys.argv[4]), int(sys.argv[5])); rho = Fr(int(sys.argv[6]), int(sys.argv[7]))
bmax = Fr(d + 1, d - 1)
gens = sorted({rs for rs in itertools.combinations_with_replacement([Fr(0), bmax, rE], n)
               if sum(1 for r in rs if r == rE and rE > bmax) <= m3})
def esym(rs, k):
    return sum(prod(c) for c in itertools.combinations(rs, k)) if k else Fr(1)
cols = [[esym(rs, k) / comb(n, k) for k in range(n + 1)] for rs in gens]
target = [rho ** k for k in range(n + 1)]

def solve(A, b):          # exact Gaussian elimination, A square (list of rows)
    m = len(A); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(m):
        p = next((r for r in range(c, m) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(m):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][m] / M[i][i] for i in range(m)]

found = None
for basis in itertools.combinations(range(len(gens)), n + 1):
    A = [[cols[j][k] for j in basis] for k in range(n + 1)]
    w = solve(A, target)
    if w is not None and all(x >= 0 for x in w):
        found = (basis, w); break
if found:
    print(f'n={n} d={d} m3={m3} rE={rE} rho={rho} (alpha={(rho-1)/(rho+1)}): FEASIBLE, exact certificate:')
    for j, wj in zip(*found):
        if wj:
            print('   ', tuple(str(x) for x in gens[j]), wj)
else:
    print(f'n={n} d={d} m3={m3} rE={rE} rho={rho} (alpha={(rho-1)/(rho+1)}): infeasible over all bases')
