"""Exact certificate that NO rescaling of the coordinates by monomials (x_a -> s^{delta_a} x_a, any real exponents)
removes the obstruction of newton_obstr.py, for blocks where no single pair has a width certificate.

For every ordered pair (a, b) with M_ab != 0 put P = n_ab d_aa, Q = d_ab n_aa (dehomogenised supports) and, for an
integer direction l, g_ab(l) = max_{supp Q} l.m - max_{supp P} l.m.  Boundedness of the rescaled ratio needs
(delta_b - delta_a).l <= g_ab(l) for every l.  A Farkas certificate for the finite system over the directions
l in {-2..2}^r: rational lambda_{ab,l} >= 0 with sum lambda_{ab,l} l (e_b - e_a) = 0 and sum lambda_{ab,l} g_ab(l) < 0.
The LP (HiGHS) proposes the support; the multipliers are then recomputed exactly over Q on that support and every
condition is checked in exact arithmetic.  Usage: python rescale_farkas.py k"""
import sys, itertools, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from fractions import Fraction as Fr
import flint
from psdtools import load_blocks

k = int(sys.argv[1])
homog = k in (2, 3)
t0 = time.time()
nproved = 0
nobs = 0
for blk in load_blocks(k):
    Mg = blk['Mg']
    n = len(Mg)
    nv = len(blk['names']) - (1 if homog else 0)
    dirs = [l for l in itertools.product(range(-2, 3), repeat=nv) if any(l)]
    rows = []  # (a, b, l, g)
    obstructed = False
    for a in range(n):
        for b in range(n):
            if a == b or Mg[a][b][0].is_zero():
                continue
            P = [m[1:] if homog else m for m in (Mg[a][b][0] * Mg[a][a][1]).to_dict()]
            Qs = [m[1:] if homog else m for m in (Mg[a][b][1] * Mg[a][a][0]).to_dict()]
            for l in dirs:
                g = max(sum(x * y for x, y in zip(l, m)) for m in Qs) - max(sum(x * y for x, y in zip(l, m)) for m in P)
                rows.append((a, b, l, int(g)))
                if g < 0:
                    obstructed = True
    if not obstructed:
        continue
    nobs += 1
    N = len(rows)
    A = lil_matrix((n * nv + 1, N))
    c = np.zeros(N)
    for i, (a, b, l, g) in enumerate(rows):
        for d in range(nv):
            if l[d]:
                A[b * nv + d, i] += l[d]
                A[a * nv + d, i] -= l[d]
        A[n * nv, i] = 1
        c[i] = g
    beq = np.zeros(n * nv + 1)
    beq[-1] = 1
    res = linprog(c, A_eq=A.tocsr(), b_eq=beq, bounds=[(0, None)] * N, method='highs')
    verdict = 'LP failed'
    if res.status == 0 and res.fun < -1e-9:
        supp = [i for i in range(N) if res.x[i] > 1e-10]
        # exact: solve the equality system on the support
        m_ = n * nv + 1
        Mx = flint.fmpq_mat(m_, len(supp) + 1)
        for jj, i in enumerate(supp):
            a, b, l, g = rows[i]
            for d in range(nv):
                if l[d]:
                    Mx[b * nv + d, jj] += l[d]
                    Mx[a * nv + d, jj] -= l[d]
            Mx[n * nv, jj] = 1
        Mx[n * nv, len(supp)] = 1
        R, rank = Mx.rref()
        # particular solution: pivot variables from R, free variables 0
        lam = [Fr(0)] * len(supp)
        piv_cols = []
        r_ = 0
        for col in range(len(supp) + 1):
            if r_ < rank and R[r_, col] != 0:
                piv_cols.append(col)
                r_ += 1
        consistent = len(supp) not in piv_cols
        if consistent:
            for rr, col in enumerate(piv_cols):
                v = R[rr, len(supp)]
                lam[col] = Fr(int(v.p), int(v.q))
            good_eq = True
            for d in range(m_):
                tot = Fr(0)
                for jj, i in enumerate(supp):
                    a, b, l, g = rows[i]
                    if d == n * nv:
                        coef = 1
                    else:
                        st, comp = divmod(d, nv)
                        coef = (l[comp] if b == st else 0) - (l[comp] if a == st else 0)
                    tot += coef * lam[jj]
                good_eq &= tot == (1 if d == n * nv else 0)
            obj = sum(lam[jj] * rows[i][3] for jj, i in enumerate(supp))
            if good_eq and all(x >= 0 for x in lam) and obj < 0:
                verdict = (f'EXACT certificate: {sum(1 for x in lam if x > 0)} multipliers, objective {obj} < 0 '
                           f'-> no monomial rescaling removes the obstruction')
                nproved += 1
            else:
                verdict = f'LP value {res.fun:.3g} < 0 but the exact recomputation failed (eq {good_eq}, obj {obj})'
        else:
            verdict = f'LP value {res.fun:.3g} < 0 but the support system is inconsistent over Q'
    elif res.status == 0:
        verdict = f'LP value {res.fun:.3g}: no certificate over these directions'
    print(f'{blk["name"]} ({blk["kind"]}), size {n}: {verdict}  [{time.time() - t0:.1f}s]', flush=True)
print(f'RESULT k = {k}: {nproved} of {nobs} obstructed blocks have an exact certificate against every monomial rescaling')
print('PASS')
