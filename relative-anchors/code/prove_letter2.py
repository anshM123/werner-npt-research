"""All-k proof for the CASE I letter-2 torus blocks of (R_inf)_j at equal weights (Theorem E3, README Section 4).
Parameters q0 = j - tau - i, t = k - j - i... (t = m - q0, m = k - tau - 2i), i, tau = |T|; modes sym / anti.
Families: q0 in {-1,0,1,2*,S}, t in {-1,0,1,2,3,S}, i in {0,1,S}, tau always symbolic (sym: tau >= 1, anti: tau >= 2);
S = symbolic with q0 >= 2 (>= 3 when t = -1; the line q0 = 2, t = -1 is separate), t >= 4, i >= 2.
For each family: (0) the level offsets met must keep the Newton window exact on the family; (1) Q, G over QQ(.)
(orbit2.SymF2); (2) sanity: exact agreement with the numeric calculus Num2 at sample points; (3) kappa =
A0 + A1 + B10 - B01 + C1 (sym) lies in ker Q, one support vector dropped; (4) symmetric elimination with greedy
smallest pivots, each pivot's numerator and denominator one-signed with nonzero constant term after shifting every
symbolic parameter to its minimum; a remaining identically-zero block is allowed (reported)."""
import sys, time, itertools
from fractions import Fraction
from math import comb
from orbit2 import SymF2, Num2, reduced2

def shifts(q0, t, i, tau, mode):
    # U empty (i = 0 and m = q0 + t = 0) forces k = tau, so tau >= 2 there
    u0 = (i == 0 and q0 is not None and t is not None and q0 + t == 0)
    return {'q0': 3 if t == -1 else 2, 't': 4, 'i': 2, 'tau': 2 if (mode == 'anti' or u0) else 1}

def univariate_sign(S):
    """S: fmpq_mpoly in one variable w >= 0: sign if no real root in [0, oo) (sympy count_roots, exact), else 0"""
    import sympy as sp
    w = sp.Symbol('w')
    P = sp.Poly(sum(sp.Rational(int(c.p), int(c.q)) * w ** e[0] for e, c in S.to_dict().items()), w)
    if P.count_roots(0, None) != 0:
        return 0
    v = P.eval(0)
    return 1 if v > 0 else (-1 if v < 0 else 0)

def sign_ok(P, be, SH):
    gens = be.Kf.ctx.gens()
    S = P.compose(*[g + SH[nm] for g, nm in zip(gens, be.syms)]) if be.syms else P
    co = list(S.coeffs())
    if not co:
        return 0
    c0 = S.to_dict().get((0,) * max(1, len(be.syms)), 0) if be.syms else co[0]
    if all(x > 0 for x in co) and c0 > 0:
        return 1
    if all(x < 0 for x in co) and c0 < 0:
        return -1
    if len(be.syms) == 1:
        return univariate_sign(S)
    return 0

def evalF(x, vals):
    if x == 0:
        return Fraction(0)
    nn = x.n(*vals) if vals else x.n(0)
    dd = x.d(*vals) if vals else x.d(0)
    r = nn / dd
    return Fraction(int(r.p), int(r.q))

def prove(q0, t, i, mode, verbose=False):
    t0 = time.time()
    tau = None
    if q0 is not None and t is not None and q0 + t < 0:
        return 'empty', 0
    be = SymF2(q0, t, i, tau)
    SH = shifts(q0, t, i, tau, mode)
    names, G, Q = reduced2(be, mode)
    if not names:
        return 'empty', 0
    # (0) window exactness on the family
    if q0 is None and min(d for d, dp in be.seen) + SH['q0'] < 2:
        return 'WINDOW-q0', 0
    if t is None and min(dp for d, dp in be.seen) + SH['t'] < 2:
        return 'WINDOW-t', 0
    mmin = (SH['q0'] if q0 is None else q0) + (SH['t'] if t is None else t)
    if (q0 is None or t is None) and mmin < 2:
        return 'M-TOO-SMALL', 0
    # (2) sanity at sample points
    syms = be.syms
    for sample in itertools.product(*[(SH[s], SH[s] + 1, SH[s] + 3) for s in syms]):
        val = dict(zip(syms, sample))
        q0v = val.get('q0', q0); tv = val.get('t', t); iv = val.get('i', i); tauv = val.get('tau', tau)
        m = q0v + tv
        k = m + 2 * iv + tauv; j = q0v + iv + tauv
        n2, G2, Q2 = reduced2(Num2(k, j, tauv, iv), mode)
        N0 = 2 ** iv * comb(m, q0v + be.sh)
        if n2 != names:
            return f'SANITY-NAMES {names} {n2} at k={k},j={j},tau={tauv},i={iv}', time.time() - t0
        for a in range(len(names)):
            for b in range(len(names)):
                if evalF(Q[a][b], sample) * N0 != Q2[a][b] or evalF(G[a][b], sample) * N0 != G2[a][b]:
                    return f'SANITY-VALUES at k={k},j={j},tau={tauv},i={iv}', time.time() - t0
    n = len(names)
    kap = {'A0': 1, 'A1': 1, 'B10': 1, 'B01': -1, 'C1': 1} if mode == 'sym' else {}
    present = [nm for nm in ['A1', 'A0', 'B10', 'B01', 'C1'] if nm in names and nm in kap]
    drop = None
    if present:
        Qk = [sum((Q[r][names.index(nm)] * kap[nm] for nm in present), be.K(0)) for r in range(n)]
        if not all(x == 0 for x in Qk):
            return 'KAPPA-NOT-KERNEL', time.time() - t0
        drop = names.index(present[0])
    act = [x for x in range(n) if x != drop]
    M = [row[:] for row in Q]
    npiv = 0
    while act:
        cand = [x for x in act if M[x][x] != 0]
        if not cand:
            if all(M[a][b] == 0 for a in act for b in act):
                return f'PSD (pivots {npiv}, zero block {len(act)}, kappa {"yes" if present else "no"})', time.time() - t0
            return 'ZERO-DIAG-NONZERO-BLOCK', time.time() - t0
        x = min(cand, key=lambda z: len(M[z][z].n) + len(M[z][z].d))
        piv = M[x][x]
        sn, sd = sign_ok(piv.n, be, SH), sign_ok(piv.d, be, SH)
        if sn == 0 or sd == 0 or sn != sd:
            return f'PIVOT-SIGN-UNPROVED at {names[x]} (num {len(piv.n)}, den {len(piv.d)} terms)', time.time() - t0
        npiv += 1
        if verbose:
            print(f'   pivot {names[x]}: num {len(piv.n)} / den {len(piv.d)} terms, positive ({time.time()-t0:.1f}s)', flush=True)
        act.remove(x)
        for r in act:
            f = M[r][x] / piv
            if f != 0:
                for s in act:
                    M[r][s] = M[r][s] - f * M[x][s]
    return f'PSD (pivots {npiv}, kappa {"yes" if present else "no"})', time.time() - t0

if __name__ == '__main__':
    modes = sys.argv[1].split(',') if len(sys.argv) > 1 else ['anti', 'sym']
    fams = []
    for mode in modes:
        for q0 in (-1, 0, 1, None):
            for t in (-1, 0, 1, 2, 3, None):
                for i in (0, 1, None):
                    fams.append((mode, q0, t, i))
        for i in (0, 1, None):
            fams.append((mode, 2, -1, i))
    allok = True; T0 = time.time()
    for mode, q0, t, i in fams:
        nsym = 1 + sum(v is None for v in (q0, t, i))
        if mode == 'sym' and nsym == 4:
            print('[sym: q0>=2, t>=4, i>=2, tau>=1] -> proved separately by prove_letter2_split.py (irrep split)', flush=True)
            continue
        res, dt = prove(q0, t, i, mode, verbose=(nsym == 4))
        ok = res.startswith('PSD') or res == 'empty'
        allok &= ok
        SH = shifts(q0, t, i, None, mode)
        lab = lambda v, nm: f'{nm}>={SH[nm]}' if v is None else f'{nm}={v}'
        print(f'[{mode}: {lab(q0,"q0")}, {lab(t,"t")}, {lab(i,"i")}, tau>={SH["tau"]}] {res}  ({dt:.1f}s)', flush=True)
    print('ALL LETTER-2 CASE I FAMILIES PROVED' if allok else 'SOME FAMILIES NOT PROVED', f'{time.time()-T0:.0f}s')
