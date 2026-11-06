"""Families symbolic in the sizes, with sigma free.  A family fixes each of q0 = j - i in {-1, 0, 1, 2, 3},
t = u - j - i in {-2, -1, 0, 1, 2, 3} and i in {0, 1}, or leaves it symbolic (then q0 = 4 + x, t = 4 + y, i = 2 + z
with x, y, z >= 0); sigma > 0 is always symbolic.  For each family with at least one symbolic size:
  (1) the reduced form over Q(sizes, sigma) (famlib.reduced);
  (2) sanity: at sample points it equals N0 times the exact numeric form of mlib (which asserts its Newton window);
  (3) the known kernel vectors are annihilated identically;
  (4) gauge: one deleted coordinate per kernel vector, with a constant nonzero determinant;
  (5) symmetric elimination on the rest, pivots chosen greedily (smallest diagonal entry that passes the test; the
      pivots of any order are the ratios of the leading minors of a symmetric permutation); every pivot, after the
      shifts, is a ratio of two
      polynomials in (x, y, z, sigma) whose coefficients all have one sign and whose part free of x, y, z is not zero,
      so it is > 0 for all x, y, z >= 0 and sigma > 0; every entry denominator passes the same test.
Usage: python prove_fam.py MINS MAXS [bscale]: the families with MINS..MAXS symbolic sizes"""
import sys, time, itertools, pickle
from fractions import Fraction as Fr
import flint
import famlib, mlib

SMIN = famlib.SMIN


def shifted(P, C):
    ctx = C.K.ctx
    gens = ctx.gens()
    args = [g + SMIN[nm] if nm in SMIN else g for g, nm in zip(gens, C.syms)]
    return P.compose(*args)


def one_signed(P, C):
    """+1 / -1 if after the shift all coefficients have that sign and the part free of the size variables is not
    zero; 0 otherwise"""
    S = shifted(P, C)
    d = S.to_dict()
    if not d:
        return 0
    co = list(d.values())
    nsz = len(C.syms) - 1
    base = [c for e, c in d.items() if all(x == 0 for x in e[:nsz])]
    if all(c > 0 for c in co) and base:
        return 1
    if all(c < 0 for c in co) and base:
        return -1
    return 0


def evalF(x, vals):
    if x == 0:
        return Fr(0)
    vals = [flint.fmpq(Fr(v).numerator, Fr(v).denominator) for v in vals]
    r = x.n(*vals) / x.d(*vals)
    return Fr(int(r.p), int(r.q))


def fam_kernel(names, C):
    s = C.sig
    out = []
    for y0 in (0, 1):
        v = {}
        for nm, cf in (('O%d' % y0, 1), ('A0%d' % y0, 1), ('A1%d' % y0, 1), ('B10%d' % y0, 1), ('B01%d' % y0, -1)):
            if nm in names:
                v[nm] = C.K(cf)
        if v:
            out.append(('p x h%d' % y0, v, 'O%d' % y0))
    ft, fi = C.ft, C.fi
    u = C.u
    if ft == 1 and fi == 0:
        out.append(('K2', {nm: cf for nm, cf in (('G0', s), ('G1', C.one)) if nm in names}, 'G1'))
    if ft == 0 and fi == 0:
        v = {nm: cf for nm, cf in (('O0', u * s), ('O1', u), ('A10', -s * s), ('A01', -s), ('A11', -s)) if nm in names}
        out.append(('T0', v, 'A10'))
    if ft == -1 and fi == 1:
        v = {nm: cf for nm, cf in (('B110', s), ('B011', C.one), ('B101', C.one), ('B111', C.one)) if nm in names}
        out.append(('T1', v, 'B111'))
    return out


def prove(q0, t, i, verbose=False, bscale=Fr(1)):
    t0 = time.time()
    if q0 is not None and t is not None and q0 + t < 0:
        return 'empty', 0, {}
    C = famlib.FamCalc(q0, t, i, bscale=bscale)
    names, Q = famlib.reduced(C)
    if not names:
        return 'empty', time.time() - t0, {}
    syms = C.syms
    sz = [nm for nm in syms if nm != 's']
    # (2) sanity at sample points
    for sample in itertools.product(*[(SMIN[nm], SMIN[nm] + 1, SMIN[nm] + 3) for nm in sz]):
        val = dict(zip(sz, sample))
        q0v, tv, iv = val.get('q0', q0), val.get('t', t), val.get('i', i)
        m = q0v + tv
        u = m + 2 * iv
        k, j = u + 1, q0v + iv
        N0 = famlib.norm0(C, q0v, tv, iv)
        for sg in (Fr(2, 7), Fr(5)):
            n2, nr2, Q2 = mlib.reduced(mlib.Calc(k, j, iv, sigma=sg, bscale=bscale))
            if n2 != names:
                return f'SANITY-NAMES {names} {n2} at {(k, j, iv)}', time.time() - t0, {}
            pt = list(sample) + [sg]
            for a in range(len(names)):
                for b in range(len(names)):
                    if evalF(Q[a][b], pt) * N0 != Fr(int(Q2[a][b].p), int(Q2[a][b].q)):
                        return f'SANITY-VALUES at {(k, j, iv, sg)} {names[a]},{names[b]}', time.time() - t0, {}
    n = len(names)
    for a in range(n):
        for b in range(n):
            if Q[a][b] != 0 and one_signed(Q[a][b].d, C) == 0:
                return f'ENTRY DENOMINATOR NOT ONE-SIGNED at {names[a]},{names[b]}', time.time() - t0, {}
    # (3) kernel
    kv = fam_kernel(names, C) if bscale == 1 else fam_kernel(names, C)
    if bscale == 1:
        for lab, v, _ in kv:
            for a in range(n):
                s = C.zero
                for nm, cf in v.items():
                    s = s + Q[a][names.index(nm)] * cf
                if s != 0:
                    return f'KAPPA {lab} NOT IN KERNEL', time.time() - t0, {}
    # (4) gauge: delete the designated coordinate of each kernel vector; the determinant must be a nonzero constant
    S = [names.index(dl) for _, _, dl in kv]
    sub = [[v.get(names[c], C.zero) for c in S] for _, v, _ in kv]
    det = det_small(sub, C)
    if det == 0 or one_signed(det.n, C) == 0 or one_signed(det.d, C) == 0:
        return f'GAUGE DETERMINANT {det} NOT SIGN-DEFINITE', time.time() - t0, {}
    keep = [c for c in range(n) if c not in S]
    M = {(a, b): Q[a][b] for a in keep for b in keep}
    npiv = 0
    info = {'names': names, 'kernel': [lab for lab, _, _ in kv], 'gauge': [names[c] for c in S], 'maxterms': 0,
            'order': []}
    act = list(keep)
    while act:
        # greedy order: the smallest diagonal entry that passes the sign test (any order of positive pivots proves
        # positive definiteness: they are the ratios of the leading minors of a symmetric permutation of the matrix)
        cands = sorted(act, key=lambda z: (len(M[z, z].n) + len(M[z, z].d)) if M[z, z] != 0 else 10 ** 12)
        chosen = None
        for x in cands:
            piv = M[x, x]
            if piv == 0:
                continue
            sn, sd = one_signed(piv.n, C), one_signed(piv.d, C)
            if sn != 0 and sd != 0 and sn == sd:
                chosen = x
                break
        if chosen is None:
            x = cands[0]
            return (f'PIVOT-SIGN-UNPROVED (no remaining diagonal entry passes; smallest {names[x]}: '
                    f'num {len(M[x, x].n)} terms, den {len(M[x, x].d)})', time.time() - t0, info)
        x = chosen
        piv = M[x, x]
        info['maxterms'] = max(info['maxterms'], len(piv.n), len(piv.d))
        info['order'].append(names[x])
        npiv += 1
        if verbose:
            print(f'    pivot {names[x]}: num {len(piv.n)} terms, den {len(piv.d)} terms ({time.time()-t0:.1f}s)',
                  flush=True)
        act.remove(x)
        for ia, a in enumerate(act):
            if M[a, x] == 0:
                continue
            f = M[a, x] / piv
            for b in act[ia:]:
                if M[x, b] != 0:
                    M[a, b] = M[a, b] - f * M[x, b]   # symmetric: compute once, mirror
                    M[b, a] = M[a, b]
    return f'PD on complement (pivots {npiv}, kernel {len(kv)})', time.time() - t0, info


def det_small(M, C):
    n = len(M)
    if n == 0:
        return C.one
    if n == 1:
        return M[0][0]
    s = C.zero
    for c in range(n):
        if M[0][c] == 0:
            continue
        sub = [row[:c] + row[c + 1:] for row in M[1:]]
        term = M[0][c] * det_small(sub, C)
        s = s + term if c % 2 == 0 else s - term
    return s


if __name__ == '__main__':
    mins = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    maxs = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    bscale = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1)
    fams = []
    for q0 in (-1, 0, 1, 2, 3, None):
        for t in (-2, -1, 0, 1, 2, 3, None):
            for i in (0, 1, None):
                ns = sum(v is None for v in (q0, t, i))
                if max(1, mins) <= ns <= maxs:
                    fams.append((ns, q0, t, i))
    # within the same number of symbolic sizes: families with i symbolic first, then those with q0 fixed
    fams.sort(key=lambda z: (z[0], z[3] is not None, z[1] is None))
    allok = True
    T0 = time.time()
    nproved = nempty = 0
    results = {}
    for ns, q0, t, i in fams:
        res, dt, info = prove(q0, t, i, verbose=(ns >= 2), bscale=bscale)
        ok = res.startswith('PD') or res == 'empty'
        allok &= ok
        nproved += res.startswith('PD')
        nempty += res == 'empty'
        results[(q0, t, i)] = res
        lab = lambda v, nm: f'{nm}>={SMIN[nm]}' if v is None else f'{nm}={v}'
        print(f'[{lab(q0, "q0")}, {lab(t, "t")}, {lab(i, "i")}] {res}  ({dt:.1f}s)'
              + (f'  largest pivot {info["maxterms"]} terms' if info.get('maxterms') else ''), flush=True)
    if bscale == 1:
        with open(f'logs/fam_{mins}_{maxs}.pkl', 'wb') as fh:
            pickle.dump(results, fh)
    print(f'{nproved} families proved, {nempty} empty, of {len(fams)}; {time.time() - T0:.0f}s')
    print('PASS' if allok else 'FAIL')
