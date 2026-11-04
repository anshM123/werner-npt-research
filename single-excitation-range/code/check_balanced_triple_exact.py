"""EXACT (rational arithmetic) check of the bath-reduced relative-anchor conditions (README, Lemma 13) for the
balanced three-copy group vector p_G = |100> + |010> + |001> (unnormalised; R and the signs do not depend on scale):
  (D)     D_{p_G} > 0;
  (R_inf) Y = A - B compressed to {row in H_1, col in H_j} + {row = vacuum, col in H_{j+1}} is PSD, j = -1..3
          (j = -1: the single state row = col = vacuum);
  (R2)    3A + B compressed to {row = vacuum, col in H_j} is PSD, j = 0..3;
with A = M^Gamma, B = R_{p_G}^Gamma on three copies.  PSD is decided exactly by symmetric Gaussian elimination over Q."""
from fractions import Fraction as Fr
import itertools
import time

t0 = time.time()
k = 3
half = Fr(1, 2)


def m1(a, b, ap, bp):
    return Fr(int(a == ap and b == bp)) - half * int(a == b and ap == bp)


def Mel(row, col, rowp, colp):
    v = Fr(1)
    for s in range(k):
        v *= m1(row[s], col[s], rowp[s], colp[s])
        if v == 0:
            return Fr(0)
    return v


def q(st):
    return sum(x != 0 for x in st)


allst = list(itertools.product(range(3), repeat=k))
p = {(1, 0, 0): Fr(1), (0, 1, 0): Fr(1), (0, 0, 1): Fr(1)}


def psd_exact(Mx):
    """exact PSD test of a symmetric rational matrix (list of lists)."""
    M = [row[:] for row in Mx]
    act = list(range(len(M)))
    while act:
        if any(M[i][i] < 0 for i in act):
            return False
        for i in act:
            if M[i][i] == 0 and any(M[i][j] != 0 for j in act):
                return False
        piv = [i for i in act if M[i][i] > 0]
        if not piv:
            return True
        i = piv[0]
        rest = [j for j in act if j != i]
        for j in rest:
            if M[j][i] == 0:
                continue
            f = M[j][i] / M[i][i]
            for l in rest:
                if M[i][l] != 0:
                    M[j][l] -= f * M[i][l]
        act = rest
    return True


def inv_exact(Mx):
    n = len(Mx)
    A = [row[:] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(Mx)]
    for c in range(n):
        r = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[r] = A[r], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


# D on the 27 column strings, block diagonal by torus weight
def wcol(y):
    return (sum(v == 1 for v in y),) + tuple(int(v == 2) for v in y)


D = {(y, yp): sum(pa * pb * Mel(a, y, b, yp) for a, pa in p.items() for b, pb in p.items()) for y in allst for yp in allst}
blocks = {}
for y in allst:
    blocks.setdefault(wcol(y), []).append(y)
assert all(D[(y, yp)] == 0 for y in allst for yp in allst if wcol(y) != wcol(yp))
Dinv = {}
D_ok = True
for wkey, ys in blocks.items():
    Bm = [[D[(y, yp)] for yp in ys] for y in ys]
    D_ok &= psd_exact(Bm) and all(True for _ in [0])
    # positive definiteness: PSD and invertible
    Bi = inv_exact(Bm)
    for i, y in enumerate(ys):
        for j, yp in enumerate(ys):
            Dinv[(y, yp)] = Bi[i][j]
print('(D) D_{p_G} > 0 (exact, every weight block PSD and invertible):', D_ok, f'[{time.time() - t0:.1f}s]', flush=True)
rows_needed = [st for st in allst if q(st) <= 1]
MV = {}
for b in rows_needed:
    for y in allst:
        for z in allst:
            v = sum(Mel(b, y, a, z) * pa for a, pa in p.items())
            if v != 0:
                MV[(b, y, z)] = v
MVrow = {}
for (b, y, z), v in MV.items():
    MVrow.setdefault((b, y), {})[z] = v


def Rel_(b, yp, a, y):
    """R[(b, y'), (a, y)] = sum_{z,z'} MV[(b,y'),z] Dinv[z,z'] MV[(a,y),z']."""
    r1 = MVrow.get((b, yp), {})
    r2 = MVrow.get((a, y), {})
    s = Fr(0)
    for z, v1 in r1.items():
        for zp, v2 in r2.items():
            dv = Dinv.get((z, zp), 0)
            if dv:
                s += v1 * dv * v2
    return s


def Aent(a, yp, b, y):      # A[(a,y'),(b,y)] = M[(b,y'),(a,y)]
    return Mel(b, yp, a, y)


def Bent(a, yp, b, y):      # B[(a,y'),(b,y)] = R[(b,y'),(a,y)]
    return Rel_(b, yp, a, y)


vac = (0, 0, 0)
H = {j: [st for st in allst if q(st) == j] for j in range(k + 1)}
ok = D_ok
for j in range(-1, k + 1):
    S = ([(a, y) for a in H[1] for y in H[j]] if j >= 0 else []) + ([(vac, y) for y in H[j + 1]] if j + 1 <= k else [])
    Y = [[Aent(a, yp, b, y) - Bent(a, yp, b, y) for (b, y) in S] for (a, yp) in S]
    sym = all(Y[i][l] == Y[l][i] for i in range(len(S)) for l in range(len(S)))
    r_inf = sym and psd_exact(Y)
    G0 = [(vac, y) for y in H[j]] if j >= 0 else []
    C = [[3 * Aent(a, yp, b, y) + Bent(a, yp, b, y) for (b, y) in G0] for (a, yp) in G0]
    r2 = psd_exact(C) if G0 else True
    ok &= r_inf and r2
    print(f'  j = {j}: (R_inf) on {len(S)} states (symmetric: {sym}): {r_inf};  (R2) on {len(G0)} states: {r2}  [{time.time() - t0:.1f}s]', flush=True)
# sanity: the exact routine must reject a false statement: (R_inf) with Y replaced by A - (11/10) B at j = 1
S = [(a, y) for a in H[1] for y in H[1]] + [(vac, y) for y in H[2]]
Yf = [[Aent(a, yp, b, y) - Fr(11, 10) * Bent(a, yp, b, y) for (b, y) in S] for (a, yp) in S]
rej = not psd_exact(Yf)
print('sanity: the routine rejects "A - (11/10) B >= 0" on the j = 1 set:', rej)
ok &= rej
print('PASS' if ok else 'FAIL')
