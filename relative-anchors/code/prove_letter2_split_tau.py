"""Generic sym letter-2 family [q0 >= 2, t >= 4, i >= 2, tau >= 1], with the S_U-irrep split to keep sizes down.
W2 = common kernel (inside span{A_x, B_xy}) of the four S_U-equivariant functionals
   Psi -> sum_{r in U} n_r psi_r,  sum_{r in U} (1 - n_r) psi_r,  sum_{r in U} sigma_r^- psi_r,  sum_{r in U} sigma_r^+ psi_r
(each lands in a two-row module, where the (H_i, chi_i)-invariants are one-dimensional); generically W2 is the
multiplicity space of (u-i-1, i, 1).  V8 = span{X1..X4, G0, C1, D0, D1} (integer basis of the (u-i, i) part, see
below).  P = [V8 basis, W2 basis].
The proof only uses: Q' = P^T Q P is block diagonal (checked identically), every denominator in P and det P are
one-signed after the shift (so P(x) is defined and invertible at every point of the family), and both blocks
eliminate with one-signed pivots (the 8-block after dropping kappa)."""
import sys, time
from orbit2 import SymF2, reduced2, Frame2, frame_labels, SYM
import prove_letter2 as PL

t0 = time.time()
q0f = None if len(sys.argv) < 2 or sys.argv[1] == 'S' else int(sys.argv[1])
tf = None if len(sys.argv) < 3 or sys.argv[2] == 'S' else int(sys.argv[2])
if_ = None if len(sys.argv) < 4 or sys.argv[3] == 'S' else int(sys.argv[3])
tauf = None if len(sys.argv) < 5 or sys.argv[4] == "S" else int(sys.argv[4])
be = SymF2(q0f, tf, if_, tauf)
SH = PL.shifts(q0f, tf, if_, None, 'sym')
names, G, Q = reduced2(be, 'sym')
assert names == SYM, names
# (0) window exactness and (2) sanity against Num2 at sample points (as in prove_letter2.prove)
import itertools
from math import comb
from orbit2 import Num2
dmin = min(d for d, dp in be.seen); dpmin = min(dp for d, dp in be.seen)
assert (q0f is not None or dmin + SH['q0'] >= 2) and (tf is not None or dpmin + SH['t'] >= 2), 'window'
for sample in itertools.product(*[(SH[s], SH[s] + 1, SH[s] + 3) for s in be.syms]):
    val = dict(zip(be.syms, sample))
    q0v = val.get('q0', q0f); tv = val.get('t', tf); iv = val.get('i', if_); tauv = val.get('tau', tauf)
    m = q0v + tv; k = m + 2 * iv + tauv; j = q0v + iv + tauv
    n2, G2, Q2 = reduced2(Num2(k, j, tauv, iv), 'sym')
    N0 = 2 ** iv * comb(m, q0v + be.sh)
    assert n2 == names
    for a in range(10):
        for b in range(10):
            assert PL.evalF(Q[a][b], sample) * N0 == Q2[a][b] and PL.evalF(G[a][b], sample) * N0 == G2[a][b], 'sanity'
print('window exact and sanity at sample points: PASS', flush=True)
K = be.K
print('built', f'{time.time()-t0:.1f}s', flush=True)
AB = ['A0', 'A1', 'B00', 'B01', 'B10', 'B11']


def sp_op(fr, lab, v):
    s = fr.pos[lab]
    return {(L[:s] + (1,) + L[s + 1:], d): c for (L, d), c in v.items() if L[s] == 0}


F = [[K(0)] * 6 for _ in range(4)]
for (s, sz) in [('R0', be.m), ('P0a', 2 * be.i)]:
    fr = Frame2(be, frame_labels(s))
    hj, hjm, hjp = fr.expand({}, 0), fr.expand({}, -1), fr.expand({}, 1)
    for c, nm in enumerate(AB):
        psi = fr.psi(nm, (s, 1))
        if not psi:
            continue
        npsi = fr.nop(s, psi)
        F[0][c] = F[0][c] + sz * fr.dot(hj, npsi)
        F[1][c] = F[1][c] + sz * (fr.dot(hj, psi) - fr.dot(hj, npsi))
        F[2][c] = F[2][c] + sz * fr.dot(hjm, fr.sm(s, 1, psi))
        F[3][c] = F[3][c] + sz * fr.dot(hjp, sp_op(fr, s, psi))


def nullspace(M, ncols):
    """reduced row echelon over the field; returns basis of the nullspace and pivot columns"""
    R = [row[:] for row in M]
    piv = []
    r = 0
    for c in range(ncols):
        pr = next((x for x in range(r, len(R)) if R[x][c] != 0), None)
        if pr is None:
            continue
        R[r], R[pr] = R[pr], R[r]
        pv = R[r][c]
        R[r] = [x / pv for x in R[r]]
        for x in range(len(R)):
            if x != r and R[x][c] != 0:
                f = R[x][c]
                R[x] = [a - f * b for a, b in zip(R[x], R[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fc in free:
        v = [K(0)] * ncols
        v[fc] = K(1)
        for row, pc in zip(R, piv):
            v[pc] = -row[fc]
        basis.append(v)
    return basis, piv


Wab, pivF = nullspace(F, 6)
print('rank of functionals', len(pivF), 'nullspace dim', len(Wab), f'{time.time()-t0:.1f}s', flush=True)
assert len(Wab) == 2
idx = {nm: names.index(nm) for nm in names}
W = []
for w in Wab:
    v = [K(0)] * 10
    for c, nm in enumerate(AB):
        v[idx[nm]] = w[c]
    W.append(v)
# the (u-i, i) part, with an integer basis: X1 = A1 + B10 = sum_r e_r (x) n_r h_j, X2 = A0 - B01 = sum_r e_r (x)
# (1 - n_r) h_j, X3 = A1 - B11 = sum_r e_r (x) sigma_r^+ h_{j-1}, X4 = A0 + B00 = sum_r e_r (x) sigma_r^- h_{j+1}
# (images of two-row twisted vectors under S_U-equivariant maps, so they lie in the (u-i, i) isotypic part),
# and the coordinate vectors G0, C1, D0, D1 (rows in T or vacuum: two-row U-part only).
def coord(**kw):
    v = [K(0)] * 10
    for nm, c in kw.items():
        v[idx[nm]] = K(c)
    return v
V8 = [coord(A1=1, B10=1), coord(A0=1, B01=-1), coord(A1=1, B11=-1), coord(A0=1, B00=1),
      coord(G0=1), coord(C1=1), coord(D0=1), coord(D1=1)]
V8names = ['X1', 'X2', 'X3', 'X4', 'G0', 'C1', 'D0', 'D1']
P = V8 + W                              # columns
# block diagonality
def bil(u, v):
    return sum((u[x] * sum((Q[x][y] * v[y] for y in range(10)), K(0)) for x in range(10)), K(0))
Qp = [[bil(P[a], P[b]) for b in range(10)] for a in range(10)]
offd = all(Qp[a][b] == 0 for a in range(8) for b in range(8, 10))
print('block diagonal (8 + 2):', offd, f'{time.time()-t0:.1f}s', flush=True)
assert offd
# validity of P: denominators and determinant one-signed after shift
dens = set()
for col in P:
    for x in col:
        if x != 0:
            dens.add(str(x.d))
okden = True
for col in P:
    for x in col:
        if x != 0 and PL.sign_ok(x.d, be, SH) == 0:
            okden = False
# det P by elimination
Mdet = [[P[c][r] for c in range(10)] for r in range(10)]
det = K(1)
A_ = [row[:] for row in Mdet]
for c in range(10):
    pr = next((x for x in range(c, 10) if A_[x][c] != 0), None)
    if pr is None:
        det = K(0); break
    if pr != c:
        A_[c], A_[pr] = A_[pr], A_[c]; det = -det
    det = det * A_[c][c]
    for x in range(c + 1, 10):
        if A_[x][c] != 0:
            f = A_[x][c] / A_[c][c]
            A_[x] = [a - f * b for a, b in zip(A_[x], A_[c])]
okdet = det != 0 and PL.sign_ok(det.n, be, SH) != 0 and PL.sign_ok(det.d, be, SH) != 0
print('P entries pole-free on the family:', okden, '; det P one-signed:', okdet, f'{time.time()-t0:.1f}s', flush=True)


def eliminate(M, drop=None, label=''):
    n = len(M)
    M = [row[:] for row in M]
    act = [x for x in range(n) if x != drop]
    npiv = 0
    while act:
        cand = [x for x in act if M[x][x] != 0]
        if not cand:
            return all(M[a][b] == 0 for a in act for b in act), npiv, len(act)
        x = min(cand, key=lambda z: len(M[z][z].n) + len(M[z][z].d))
        piv = M[x][x]
        sn, sd = PL.sign_ok(piv.n, be, SH), PL.sign_ok(piv.d, be, SH)
        if sn == 0 or sd == 0 or sn != sd:
            print(f'  {label}: pivot sign UNPROVED ({len(piv.n)}/{len(piv.d)} terms)')
            return False, npiv, len(act)
        npiv += 1
        print(f'  {label} pivot {npiv}: {len(piv.n)}/{len(piv.d)} terms, positive ({time.time()-t0:.1f}s)', flush=True)
        act.remove(x)
        for r in act:
            f = M[r][x] / piv
            if f != 0:
                for s in act:
                    M[r][s] = M[r][s] - f * M[x][s]
    return True, npiv, 0


Q2 = [[Qp[a][b] for b in range(8, 10)] for a in range(8, 10)]
ok2, n2, z2 = eliminate(Q2, None, '(u-i-1,i,1)')
# kappa = X1 + X2 + C1 in V8 coordinates
kap = {'A0': 1, 'A1': 1, 'B10': 1, 'B01': -1, 'C1': 1}
kv = [1, 1, 0, 0, 0, 1, 0, 0]
free = list(range(8))
Q8 = [[Qp[a][b] for b in range(8)] for a in range(8)]
Qk = [sum((Q8[a][b] * kv[b] for b in range(8)), K(0)) for a in range(8)]
okk = all(x == 0 for x in Qk)
# kappa in V8 must reproduce kappa: check P[:8] kv == kappa
recon = [sum((P[b][x] * kv[b] for b in range(8)), K(0)) for x in range(10)]
okk2 = all(recon[x] == kap.get(names[x], 0) for x in range(10))
drop = next(b for b in range(8) if kv[b] != 0)
print('kappa in ker Q8:', okk, '; kappa in V8:', okk2, '; drop', V8names[drop], flush=True)
ok8, n8, z8 = eliminate(Q8, drop, '(u-i,i)')
allok = okden and okdet and okk and okk2 and ok2 and ok8 and z2 == 0 and z8 == 0
print('GENERIC SYM FAMILY PROVED' if allok else 'NOT PROVED', f'(pivots {n2} + {n8}) {time.time()-t0:.0f}s')
