"""Reusable tests on a symmetric matrix of rational functions, given as a list of lists of (num, den) fmpq_mpoly.
  obstruction(Mg)       exact: ordered pairs with M_ab / M_aa unbounded (integer direction certificate), and how many
                        of them also have a width certificate (no monomial rescaling helps);
  mono_psd(Mg)          exact: number of monomials of den * Mg whose coefficient matrix is not PSD (q = 1);
  rescale_lp(Mg, homog) floating point LP for a monomial rescaling removing the obstruction; returns rounded exponents
                        (of s -> s^2, i.e. half-integer exponents in s) or None;
  apply_rescale(...)    the rescaled matrix in the variables u with s = u^2 (exact)."""
import itertools
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from psdtools import coeff_mats, psd_exact, lcm_poly


def poly_of(Mg):
    n = len(Mg)
    den = lcm_poly([Mg[a][c][1] for a in range(n) for c in range(n)])
    if den.leading_coefficient() < 0:
        den = -den
    return den, [[Mg[a][c][0] * (den / Mg[a][c][1]) for c in range(n)] for a in range(n)]


def obstruction(Mg, rng=3):
    n = len(Mg)
    nv = Mg[0][0][0].context().nvars()
    dirs = [l for l in itertools.product(range(-rng, rng + 1), repeat=nv) if any(l)]
    dirs.sort(key=lambda l: sum(abs(x) for x in l))
    certs, wcerts = [], []
    for a in range(n):
        na, da = Mg[a][a]
        if not all(c > 0 for c in na.coeffs()):
            return None  # diagonal not subtraction-free in reduced form; test not applicable
        for b in range(n):
            if b == a or Mg[a][b][0].is_zero():
                continue
            nab, dab = Mg[a][b]
            P = nab * da
            Q = dab * na
            if not all(c > 0 for c in Q.coeffs()):
                Q = -Q
            if not all(c > 0 for c in Q.coeffs()):
                return None
            sP, sQ = list(P.to_dict()), list(Q.to_dict())
            for l in dirs:
                vP = [sum(x * y for x, y in zip(l, m)) for m in sP]
                mp = max(vP)
                if mp > max(sum(x * y for x, y in zip(l, m)) for m in sQ) and vP.count(mp) == 1:
                    certs.append((a, b, l))
                    break
            else:
                continue
            for l in dirs:
                vP = [sum(x * y for x, y in zip(l, m)) for m in sP]
                vQ = [sum(x * y for x, y in zip(l, m)) for m in sQ]
                if max(vP) - min(vP) > max(vQ) - min(vQ):
                    wcerts.append((a, b, l))
                    break
    return certs, wcerts


def mono_psd(Mg):
    den, Nt = poly_of(Mg)
    C = coeff_mats(Nt)
    return len(C), sum(1 for m in C if not psd_exact(C[m])[0])


def rescale_lp(Mg, homog):
    n = len(Mg)
    nvall = Mg[0][0][0].context().nvars()
    nv = nvall - (1 if homog else 0)

    def dh(p):
        return [m[1:] if homog else m for m in p.to_dict()]
    pairs = []
    for a in range(n):
        for b in range(n):
            if a != b and not Mg[a][b][0].is_zero():
                pairs.append((a, b, dh(Mg[a][b][0] * Mg[a][a][1]), dh(Mg[a][b][1] * Mg[a][a][0])))
    if not pairs:
        return np.zeros((n, nv), dtype=int)
    nd = n * nv
    blocks = []
    nl = 0
    for (a, b, sP, sQ) in pairs:
        for m in sP:
            blocks.append((a, b, m, sQ, nl))
            nl += len(sQ)
    if nl > 60000:
        return 'skip'  # LP too large for the memory budget
    A = lil_matrix((len(blocks) * (nv + 1) + nv, nd + nl))
    rhs = []
    r = 0
    for (a, b, m, sQ, off) in blocks:
        for i in range(nv):
            for jq, q in enumerate(sQ):
                if q[i]:
                    A[r, nd + off + jq] = q[i]
            A[r, b * nv + i] += -1
            A[r, a * nv + i] += 1
            rhs.append(m[i])
            r += 1
        for jq in range(len(sQ)):
            A[r, nd + off + jq] = 1
        rhs.append(1)
        r += 1
    for i in range(nv):
        A[r, i] = 1
        rhs.append(0)
        r += 1
    res = linprog(np.zeros(nd + nl), A_eq=A.tocsr(), b_eq=np.array(rhs, dtype=float),
                  bounds=[(-20, 20)] * nd + [(0, None)] * nl, method='highs')
    if res.status != 0:
        return None
    d2 = np.round(2 * res.x[:nd].reshape(n, nv)).astype(int)
    return d2 - d2.min(axis=0)


def apply_rescale(Mg, d2, homog):
    """entries M_ab(u^2) u^{d2_a + d2_b} (variables u, s = u^2)."""
    n = len(Mg)
    ctx = Mg[0][0][0].context()

    def sq(p):
        return ctx.from_dict({tuple(2 * e for e in m): c for m, c in p.to_dict().items()})

    def mono(ex):
        full = ((0,) + tuple(int(x) for x in ex)) if homog else tuple(int(x) for x in ex)
        return ctx.from_dict({full: 1})
    return [[(sq(Mg[a][b][0]) * mono(d2[a] + d2[b]), sq(Mg[a][b][1])) for b in range(n)] for a in range(n)]


def farkas(Mg, homog, rng=2):
    """exact Farkas certificate that no monomial rescaling removes the obstruction (as in rescale_farkas.py).
    returns ('certificate', objective) or ('none', lp_value) or ('not obstructed', None)."""
    from fractions import Fraction as Fr
    import flint
    n = len(Mg)
    nv = Mg[0][0][0].context().nvars() - (1 if homog else 0)
    dirs = [l for l in itertools.product(range(-rng, rng + 1), repeat=nv) if any(l)]
    rows = []
    obstructed = False
    for a in range(n):
        for b in range(n):
            if a == b or Mg[a][b][0].is_zero():
                continue
            P = [m[1:] if homog else m for m in (Mg[a][b][0] * Mg[a][a][1]).to_dict()]
            Qs = [m[1:] if homog else m for m in (Mg[a][b][1] * Mg[a][a][0]).to_dict()]
            for l in dirs:
                g = int(max(sum(x * y for x, y in zip(l, m)) for m in Qs) - max(sum(x * y for x, y in zip(l, m)) for m in P))
                rows.append((a, b, l, g))
                obstructed |= g < 0
    if not obstructed:
        return ('not obstructed', None)
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
    if res.status != 0 or res.fun >= -1e-9:
        return ('none', None if res.status != 0 else res.fun)
    supp = [i for i in range(N) if res.x[i] > 1e-10]
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
    piv = []
    r_ = 0
    for col in range(len(supp) + 1):
        if r_ < rank and R[r_, col] != 0:
            piv.append(col)
            r_ += 1
    if len(supp) in piv:
        return ('none', res.fun)
    lam = [Fr(0)] * len(supp)
    for rr, col in enumerate(piv):
        v = R[rr, len(supp)]
        lam[col] = Fr(int(v.p), int(v.q))
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
        if tot != (1 if d == n * nv else 0):
            return ('none', res.fun)
    obj = sum(lam[jj] * rows[i][3] for jj, i in enumerate(supp))
    if all(x >= 0 for x in lam) and obj < 0:
        return ('certificate', obj)
    return ('none', res.fun)
