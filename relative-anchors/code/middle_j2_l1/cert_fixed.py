"""Certificate for the letter-1 block of (R_inf)_2 at fixed level sizes, for every q > 0.

python cert_fixed.py U0 U1 [bscale] [save=file]   (save: the certified forms S, as coefficient lists, for check_reduced.py)
  U0 copies of squared weight 1 and U1 copies of squared weight q (k = U0 + U1 >= 5).
Steps, all exact over Q(q) (ufield.RF):
 1. the entries (l2block.Block over Q(q));
 2. the kernel identities Y_PP 1_X = 0 and Y_QP 1_X = 0 for every pair X, identically in q (b = 1);
 3. for every irrep of S_U0 x S_U1 that occurs: the test vectors (young.py), the dimension counts
    sum dim(lam) mult(lam, M_P) = k C(k,2), sum dim(lam) mult(lam, K) = C(k,2), sum dim(lam) mult(lam, M_Q) = C(k,3),
    the reduced kernel identities, and the gauge-fixed Schur complement S of the vacuum states (fixedred.py);
 4. every pivot of the symmetric elimination of S (= ratio of consecutive leading principal minors), hence every
    leading principal minor, is > 0 for every q > 0: numerator and denominator are polynomials in q whose coefficients
    are all >= 0, or (if not) have no root in (0, inf) by an exact Sturm count and are > 0 at q = 1.
Then the block is >= 0 with kernel exactly span{p (x) |1@P> : P a pair} for every q > 0 (NOTES.md, Section 3).
With bscale != 1 (negative control: every term with E^{-1} multiplied by bscale) step 2 is skipped and the run prints
FAIL as soon as one pivot is not positive on (0, inf) (with a rational point where it is negative), PASS otherwise.
Prints PASS or FAIL on the last line."""
import sys, time, pickle
from fractions import Fraction as Fr
from math import comb
from ufield import RF, Q, rf_positive
from fixedred import build_block, reduce_all, leading_pivots

t0 = time.time()


def log(*a):
    print(*a, flush=True)


def find_negative_point(f):
    """a rational q > 0 at which the rational function f is <= 0 or undefined, or None.  Candidates: a grid, and
    rational points between and around the positive real roots of the numerator and the denominator (located in
    floating point only to choose candidates; the sign at each candidate is evaluated exactly)."""
    import numpy as np
    cands = [Fr(1, 10 ** 6), Fr(1, 10 ** 4), Fr(1, 1000), Fr(1, 100), Fr(1, 30), Fr(1, 10), Fr(1, 3), Fr(1), Fr(3),
             Fr(10), Fr(30), Fr(100), Fr(1000), Fr(10 ** 4), Fr(10 ** 6)]
    roots = []
    for p in (f.n, f.d):
        cs = [float(int(c.p)) / float(int(c.q)) for c in p.coeffs()]
        if len(cs) > 1:
            for r in np.roots(cs[::-1]):
                if r.real > 0 and abs(r.imag) <= 1e-7 * abs(r):
                    roots.append(float(r.real))
    roots = sorted(roots)
    for i, r in enumerate(roots):
        cands += [Fr(r * 0.999).limit_denominator(10 ** 9), Fr(r * 1.001).limit_denominator(10 ** 9)]
        if i + 1 < len(roots):
            cands.append(Fr((r + roots[i + 1]) / 2).limit_denominator(10 ** 12))
    for x in cands:
        if x <= 0:
            continue
        try:
            if f(x) <= 0:
                return x
        except ZeroDivisionError:
            return x
    return None


def main():
    args = [a for a in sys.argv[1:] if '=' not in a]
    opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    u0, u1 = int(args[0]), int(args[1])
    b = Fr(args[2]) if len(args) > 2 else Fr(1)
    k = u0 + u1
    one = RF(1)
    q = Q if u1 > 0 else RF(1)
    log(f'letter-1 block of (R_inf)_2, k = {k}, sizes (u0, u1) = ({u0}, {u1}), squared weights 1 and q, bscale {b}')
    blk = build_block(u0, u1, q, one, b)
    log(f'  states: {len(blk.P)} P + {len(blk.Qs)} Q; inverses of Ehat_2, Ehat_3 over Q(q) [{time.time() - t0:.1f}s]')
    parts, kern_ok = reduce_all(blk, u0, u1, b, log=log)
    log(f'  reduced forms built [{time.time() - t0:.1f}s]')
    ok = True
    if b == 1:
        log(f'  kernel identities Y_PP 1_X = 0, Y_QP 1_X = 0 for all {comb(k, 2)} pairs X, identically in q: {kern_ok}')
        ok &= kern_ok
    sP = sum(r['dim'] * r['mP'] for r in parts)
    sK = sum(r['dim'] * r['d'] for r in parts)
    sQ = sum(r['dim'] * r['nQ'] for r in parts)
    dims_ok = (sP == k * comb(k, 2) and sK == comb(k, 2) and sQ == comb(k, 3)
               and all(r['mP'] == r['d'] + r['nW'] for r in parts))
    log(f'  dimension counts: P {sP} = {k * comb(k, 2)}, K {sK} = {comb(k, 2)}, Q {sQ} = {comb(k, 3)}, '
        f'mult = kernel + gauge-fixed in every part: {dims_ok}')
    ok &= dims_ok
    saved = []
    failed_point = None
    for r in parts:
        lam = r['lam']
        if b == 1 and 'ker_red_ok' in r:
            ok &= r['ker_red_ok']
        S = r['S']
        n = len(S)
        sym = all(S[i][j] == S[j][i] for i in range(n) for j in range(n))
        ok &= sym
        piv = leading_pivots(S, one)
        good_part = len(piv) == n and all(p is not None for p in piv)
        descs = []
        minor = one
        for m, p in enumerate(piv):
            if p is None:
                good_part = False
                descs.append(f'pivot {m + 1} vanishes identically')
                break
            minor = minor * p
            g, how = rf_positive(minor)
            descs.append(f'minor {m + 1}: num deg {minor.n.degree()}, den deg {minor.d.degree()}: {how}')
            if not g:
                good_part = False
                if failed_point is None:
                    x = find_negative_point(minor)
                    failed_point = (lam, m + 1, x)
                break
        ok &= good_part
        log(f'  part {lam} (dim {r["dim"]}): gauge-fixed S of size {n}, symmetric {sym}, kernel (reduced) '
            f'{r.get("ker_red_ok", "-")}; every leading minor > 0 on (0, inf): {good_part} [{time.time() - t0:.1f}s]')
        for dsc in descs:
            log('     ' + dsc)
        saved.append((lam, r['dim'], [[([str(c) for c in x.n.coeffs()], [str(c) for c in x.d.coeffs()]) for x in row]
                                      for row in S]))
    if failed_point is not None:
        lam, m, x = failed_point
        log(f'  first failing minor: part {lam}, minor {m}; nonpositive at q = {x}')
    if 'save' in opts:
        with open(opts['save'], 'wb') as fh:
            pickle.dump(dict(u0=u0, u1=u1, b=str(b), parts=saved), fh)
    log(f'time {time.time() - t0:.1f}s')
    log('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
