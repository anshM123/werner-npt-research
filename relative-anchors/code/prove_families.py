"""All-k proof of the equal-weight letter-1 block of (R_inf)_j, family by family.
Every (k, j, i) with vectors has q0 = j - i >= -1, t = k - j - i >= -1, i >= 0.  Families:
q0 in {-1, 0, 1, S}, t in {-1, 0, 1, 2, 3, S}, i in {0, 1, S}, S = symbolic with q0 >= 2, t >= 4, i >= 2.
For each family with at least one symbolic parameter:
  (1) Q, G from the orbit calculus with the exact-node backend SymF;
  (2) sanity: exact agreement with the numeric orbit calculus at sample points of the family;
  (3) kappa = A0 + A1 + B10 - B01 (p (x) h_{j,i}) lies in ker Q; drop one of its support vectors;
  (4) symmetric elimination, greedy smallest pivot; every pivot, after q0 = 2 + x, t = 4 + y, i = 2 + z,
      has numerator and denominator with all coefficients of one common sign and nonzero constant terms
      (=> pivot > 0 on the whole family); a final all-zero block is allowed (extra kernel).
Families with no symbolic parameter have k <= 6 and are covered by num_all.py.
Usage: python prove_families.py [max_symbolic_params]"""
import sys, time, itertools
from fractions import Fraction
from math import comb
import flint
import orbit3
from symf import SymF
maxs = int(sys.argv[1]) if len(sys.argv) > 1 else 3
SHIFT = {'q0': 2, 't': 4, 'i': 2}

def shifted(P, be):
    """P: fmpq_mpoly in be.Kf.ctx; returns shifted poly."""
    ctx = be.Kf.ctx
    gens = ctx.gens()
    args = [g + SHIFT[nm] for g, nm in zip(gens, be.syms)]
    return P.compose(*args)

def sign_ok(P, be):
    S = shifted(P, be)
    co = list(S.coeffs())
    if not co:
        return 0
    c0 = S.constant_coefficient() if hasattr(S, 'constant_coefficient') else None
    if c0 is None:
        d = S.to_dict(); c0 = d.get((0,) * len(be.syms), 0)
    if all(x > 0 for x in co) and c0 > 0:
        return 1
    if all(x < 0 for x in co) and c0 < 0:
        return -1
    return 0

def evalF(x, vals):
    if x == 0:
        return Fraction(0)
    n = x.n(*vals); d = x.d(*vals)
    r = n / d
    return Fraction(int(r.p), int(r.q))

def prove(q0, t, i, verbose=False):
    global SHIFT
    # with t = -1 and q0 symbolic, m = q0 - 1; start at q0 = 3 so that m >= 2 (q0 = 2 is its own family)
    SHIFT = {'q0': 3 if t == -1 else 2, 't': 4, 'i': 2}
    t0 = time.time()
    be = SymF(q0, t, i)
    if q0 is not None and t is not None and q0 + t < 0:
        return 'empty', 0
    names, G, Q = orbit3.reduced_orbit(be)
    if not names:
        return 'empty', 0
    # sanity at sample points
    syms = be.syms
    for sample in itertools.product(*[(SHIFT[s], SHIFT[s] + 1, SHIFT[s] + 3) for s in syms]):
        val = dict(zip(syms, sample))
        q0v = val.get('q0', q0); tv = val.get('t', t); iv = val.get('i', i)
        k, j = q0v + tv + 2 * iv, q0v + iv
        m = q0v + tv
        n2, G2, Q2 = orbit3.reduced_orbit(orbit3.Num(k, j, iv))
        N0 = 2 ** iv * comb(m, q0v + be.sh)
        idx = [names.index(nm) for nm in n2]
        if len(idx) != len(n2) or any(nm not in n2 for nm in names):
            return f'SANITY-NAMES {names} {n2} at {(k,j,iv)}', time.time() - t0
        for a, x in enumerate(idx):
            for b, y in enumerate(idx):
                if evalF(Q[x][y], sample) * N0 != Q2[a][b] or evalF(G[x][y], sample) * N0 != G2[a][b]:
                    return f'SANITY-VALUES at {(k,j,iv)}', time.time() - t0
    n = len(names)
    kap = {'A0': 1, 'A1': 1, 'B10': 1, 'B01': -1}
    present = [nm for nm in kap if nm in names]
    drop = None
    if present:
        Qk = [sum((Q[r][names.index(nm)] * kap[nm] for nm in present), be.K(0)) for r in range(n)]
        if not all(x == 0 for x in Qk):
            return 'KAPPA-NOT-KERNEL', time.time() - t0
        drop = names.index(present[[ 'A1', 'A0', 'B10', 'B01'].index(present[0])] if False else ([nm for nm in ['A1','A0','B10','B01'] if nm in present][0]))
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
        sn, sd = sign_ok(piv.n, be), sign_ok(piv.d, be)
        if sn == 0 or sd == 0 or sn != sd:
            return f'PIVOT-SIGN-UNPROVED at {names[x]} (num {len(piv.n)} terms, den {len(piv.d)})', time.time() - t0
        npiv += 1
        if verbose: print(f'  pivot {names[x]}: num {len(piv.n)} terms, den {len(piv.d)} terms, positive  ({time.time()-t0:.1f}s)', flush=True)
        act.remove(x)
        for r in act:
            f = M[r][x] / piv
            if f != 0:
                for s in act:
                    M[r][s] = M[r][s] - f * M[x][s]
    return f'PSD (pivots {npiv}, kappa {"yes" if present else "no"})', time.time() - t0

fams = [(1, 2, -1, None)]   # the line q0 = 2, t = -1 (m = 1), i >= 2
for q0 in (-1, 0, 1, None):
    for t in (-1, 0, 1, 2, 3, None):
        for i in (0, 1, None):
            ns = sum(v is None for v in (q0, t, i))
            if 1 <= ns <= maxs:
                fams.append((ns, q0, t, i))
fams.sort(key=lambda z: z[0])
allok = True
T0 = time.time()
for ns, q0, t, i in fams:
    res, dt = prove(q0, t, i, verbose=(ns == 3))
    ok = res.startswith('PSD') or res == 'empty'
    allok &= ok
    lab = lambda v, nm: f'{nm}>={SHIFT[nm]}' if v is None else f'{nm}={v}'  # SHIFT as set by the last prove()
    print(f'[{lab(q0,"q0")}, {lab(t,"t")}, {lab(i,"i")}] {res}  ({dt:.1f}s)', flush=True)
print('ALL FAMILIES PROVED' if allok else 'SOME FAMILIES NOT PROVED', f'{time.time()-T0:.0f}s')
