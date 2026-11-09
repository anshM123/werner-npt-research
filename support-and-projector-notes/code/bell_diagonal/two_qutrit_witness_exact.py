"""Exact version of the two-qutrit Bell-diagonal witness (Section 7, "Bell-diagonal operators").

Construction (following an independent audit, re-implemented here):
  * weights h'_z = (round(10^6 h_z) + 100) / 10^6, rational, from data/bell_k2_witness_h.txt;
  * each entry of the floating matrix C (data/bell_k2_witness_C_*.txt) is replaced by the Eisenstein integer u + v w,
    w = e^{2 pi i/3}, u = round(10^6 (Re C + Im C/sqrt 3)), v = round(2 10^6 Im C/sqrt 3);
  * with E that integer matrix, A its first two columns, R its first two rows and M = E[0:2, 0:2], put
    C~ = A adj(M) R.  Its upper-left 2 x 2 minor is det(M)^3, so C~ has rank exactly two when det M != 0.
Everything below is exact (integers and Fractions; Eisenstein integers as pairs a + b w with w^2 = -1 - w).
It checks: rank two, full local support of range and co-range at both copies (nonzero 3 x 3 minors of flattenings),
all 1440 stabilizer-plane cone conditions for h' with a strictly positive margin, and
<vec C~ | K_{h'} | vec C~> / ||C~||^2 < -1/300.
usage: python two_qutrit_witness_exact.py"""
import itertools, os, math
from fractions import Fraction as Fr
import numpy as np

here = os.path.dirname(os.path.abspath(__file__)); data = os.path.join(here, '..', '..', 'data')
h = np.loadtxt(os.path.join(data, 'bell_k2_witness_h.txt'))
Cf = np.loadtxt(os.path.join(data, 'bell_k2_witness_C_real.txt')) + 1j * np.loadtxt(os.path.join(data, 'bell_k2_witness_C_imag.txt'))
SC = 10 ** 6
hq = [Fr(int(round(SC * x)) + 100, SC) for x in h]

# Eisenstein integers a + b w
def emul(x, y):
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)
def eadd(x, y): return (x[0] + y[0], x[1] + y[1])
def esub(x, y): return (x[0] - y[0], x[1] - y[1])
def econj(x): return (x[0] - x[1], -x[1])
def enorm(x): return x[0] * x[0] - x[0] * x[1] + x[1] * x[1]       # |a + b w|^2
def wpow(k):
    k %= 3
    return [(1, 0), (0, 1), (-1, -1)][k]
ZERO = (0, 0)

E = [[(int(round(SC * (z.real + z.imag / math.sqrt(3)))), int(round(2 * SC * z.imag / math.sqrt(3)))) for z in row] for row in Cf]
A = [[E[i][j] for j in range(2)] for i in range(9)]
R = [[E[i][j] for j in range(9)] for i in range(2)]
M = [[E[i][j] for j in range(2)] for i in range(2)]
detM = esub(emul(M[0][0], M[1][1]), emul(M[0][1], M[1][0]))
assert detM != ZERO
adj = [[M[1][1], (-M[0][1][0], -M[0][1][1])], [(-M[1][0][0], -M[1][0][1]), M[0][0]]]
AadjM = [[eadd(emul(A[i][0], adj[0][k]), emul(A[i][1], adj[1][k])) for k in range(2)] for i in range(9)]
Ct = [[eadd(emul(AadjM[i][0], R[0][j]), emul(AadjM[i][1], R[1][j])) for j in range(9)] for i in range(9)]
minor = esub(emul(Ct[0][0], Ct[1][1]), emul(Ct[0][1], Ct[1][0]))
d3 = emul(emul(detM, detM), detM)
assert minor == d3, "upper-left minor should be det(M)^3"
print(f"det M = {detM[0]} + {detM[1]} w (nonzero); upper-left 2x2 minor of C~ = det(M)^3: rank of C~ is exactly 2")

def det3(m):
    t = ZERO
    for (p, sgn) in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
        term = emul(emul(m[0][p[0]], m[1][p[1]]), m[2][p[2]])
        t = eadd(t, term) if sgn == 1 else esub(t, term)
    return t
def local_support_ok(vecs):
    """vecs: list of 9-vectors (Eisenstein); check full local support at both copies: the 3 x (3*len) flattenings
    have a nonzero 3 x 3 minor."""
    ok = []
    for site in range(2):
        cols = []
        for v in vecs:
            T = [[v[3 * a + b] for b in range(3)] for a in range(3)]          # T[a][b], a = copy 0, b = copy 1
            if site == 0: cols += [[T[a][b] for a in range(3)] for b in range(3)]
            else: cols += [[T[a][b] for b in range(3)] for a in range(3)]
        found = False
        for c3 in itertools.combinations(range(len(cols)), 3):
            m = [[cols[c][r] for c in c3] for r in range(3)]
            if det3(m) != ZERO: found = True; break
        ok.append(found)
    return ok
range_vecs = [[AadjM[i][k] for i in range(9)] for k in range(2)]
corange_vecs = [[econj(R[k][j]) for j in range(9)] for k in range(2)]
print("full local support of range at copies 0, 1:", local_support_ok(range_vecs),
      "; of co-range:", local_support_ok(corange_vecs))

# Weyl coefficients Tr(D_z^dagger C~), D_z = X^a1 Z^b1 (x) X^a2 Z^b2, (X^a Z^b)_{ij} = [i = j + a] w^{b j}
pts = list(itertools.product(range(3), repeat=4))        # (a1, a2, b1, b2)
norms = []
for (a1, a2, b1, b2) in pts:
    c = ZERO
    for j1 in range(3):
        for j2 in range(3):
            i1, i2 = (j1 + a1) % 3, (j2 + a2) % 3
            # conj(D_z)_{(i1 i2),(j1 j2)} = w^{-(b1 j1 + b2 j2)}
            c = eadd(c, emul(wpow(-(b1 * j1 + b2 * j2)), Ct[3 * i1 + i2][3 * j1 + j2]))
    norms.append(enorm(c))
total_norm = sum(enorm(Ct[i][j]) for i in range(9) for j in range(9))
assert sum(norms) == 9 * total_norm, "Weyl Parseval"
value = sum(hq[k] * norms[k] for k in range(81)) / (9 * total_norm)
print(f"Weyl Parseval exact; all 81 coefficients nonzero: {all(n_ > 0 for n_ in norms)}")
print(f"<vec C~|K_h'|vec C~>/||C~||^2 = {float(value):.13f}  (exact rational; < -1/300: {value < Fr(-1, 300)})")

# 1440 stabilizer-plane cone conditions for h' (exact)
key = {p: i for i, p in enumerate(pts)}
def symp(u, v): return (u[0] * v[2] + u[1] * v[3] - u[2] * v[0] - u[3] * v[1]) % 3
nz = [p for p in pts if any(p)]
lagr = []
for u in nz:
    for v in nz:
        if symp(u, v): continue
        span = frozenset(tuple((i * a + j * b) % 3 for a, b in zip(u, v)) for i in range(3) for j in range(3))
        if len(span) == 9 and span not in lagr: lagr.append(span)
assert len(lagr) == 40
count = 0; min_e2 = None; min_tot = None
for Ls in lagr:
    subs = []
    for w_ in Ls:
        if not any(w_): continue
        sp = frozenset(tuple((c * x) % 3 for x in w_) for c in range(3))
        if sp not in subs: subs.append(sp)
    seen = set()
    for x in pts:
        Lx = frozenset(tuple((a + b) % 3 for a, b in zip(z, x)) for z in Ls)
        if Lx in seen: continue
        seen.add(Lx)
        for I in subs:
            rem = set(Lx); sums = []
            while rem:
                z0 = next(iter(rem)); g = {tuple((a + b) % 3 for a, b in zip(z0, i_)) for i_ in I}
                sums.append(sum(hq[key[z]] for z in g)); rem -= g
            tot = sum(sums); e2 = sums[0] * sums[1] + sums[1] * sums[2] + sums[0] * sums[2]
            count += 1
            min_e2 = e2 if min_e2 is None else min(min_e2, e2)
            min_tot = tot if min_tot is None else min(min_tot, tot)
print(f"{count} cone conditions: smallest total {float(min_tot):.6f} > 0, smallest e_2 = {float(min_e2):.3e} > 0: "
      f"{min_tot > 0 and min_e2 > 0}")
