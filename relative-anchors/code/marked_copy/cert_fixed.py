"""Certificate, for each k in [KMIN, KMAX], every j = -1..k and every i = 0..floor((k-1)/2), that the reduced form Q_i
(mlib.reduced, symbolic in sigma) is positive semidefinite for every sigma > 0 with kernel exactly the span of the
known kernel vectors (mlib.kernel_vectors).  Steps per (k, j, i):
  1. Q_i over Q(sigma); every entry denominator > 0 on (0, inf); Q_i symmetric;
  2. Q_i kappa = 0 identically for every known kernel vector kappa;
  3. gauge: a set S of coordinates, one per kernel vector, with det(kappa_l[S]) a nonzero polynomial with no root in
     (0, inf) (sign test); the complement W of S is then a complement of span(kappa) for every sigma > 0;
  4. symmetric elimination of Q_i restricted to W in the order of mlib.NAMES; every pivot (a reduced rational function
     of sigma) is > 0 on (0, inf) (signtest.positive_on_halfline: nonnegative coefficients, else a witness search,
     else Sturm).  Then every leading principal minor is > 0 there (Sylvester).
Usage: python cert_fixed.py KMIN KMAX [bscale] [save.pkl]
With bscale != 1 (negative control) the kernel identities are not required; the run reports FAIL with a witness
sigma where a pivot of the gauge-fixed form is <= 0 (expected for 1001/1000)."""
import sys, time, itertools, pickle
from fractions import Fraction as Fr
import mlib
from ufield import RF
from signtest import positive_on_halfline
import flint

PL = flint.fmpq_poly


def poly_sign_definite(p):
    """+1 if p > 0 on (0, inf), -1 if p < 0 there, 0 otherwise (exact, via the sign test)"""
    f = RF(p)
    ok, how, wit = positive_on_halfline(f)
    if ok:
        return 1
    ok, how, wit = positive_on_halfline(RF(-p))
    if ok:
        return -1
    return 0


def det_poly(M):
    """determinant of a small square matrix of RF (Laplace), as RF"""
    n = len(M)
    if n == 0:
        return RF(1)
    if n == 1:
        return M[0][0]
    s = RF(0)
    for c in range(n):
        if M[0][c].is_zero():
            continue
        sub = [row[:c] + row[c + 1:] for row in M[1:]]
        term = M[0][c] * det_poly(sub)
        s = s + term if c % 2 == 0 else s - term
    return s


def choose_gauge(names, kvecs):
    """columns S (indices into names), one per kernel vector, with det sign-definite on (0, inf)"""
    r = len(kvecs)
    if r == 0:
        return [], 'none'
    rows = [[RF(v.get(nm, 0)) if not isinstance(v.get(nm, 0), RF) else v.get(nm, 0) for nm in names] for _, v in kvecs]
    # preferred order of columns: the kernel-vector names first in the order O, A, B
    pref = sorted(range(len(names)), key=lambda c: (0 if names[c][0] == 'O' else 1 if names[c][0] == 'A' else 2
                                                     if names[c][0] == 'B' else 3, c))
    for S in itertools.combinations(pref, r):
        sub = [[row[c] for c in S] for row in rows]
        d = det_poly(sub)
        if d.is_zero():
            continue
        if poly_sign_definite(d.n) != 0 and poly_sign_definite(d.d) != 0:
            return list(S), f'det {d}'
    return None, 'no gauge'


def certify(k, j, i, bscale=1, verbose=False):
    C = mlib.Calc(k, j, i, sigma=None, bscale=bscale)
    names, norms, Q = mlib.reduced(C)
    info = {'k': k, 'j': j, 'i': i, 'names': names}
    if not names:
        info['status'] = 'empty'
        return True, info
    n = len(names)
    for a in range(n):
        for b in range(n):
            if Q[a][b] != Q[b][a]:
                info['status'] = 'NOT SYMMETRIC'
                return False, info
            if not Q[a][b].is_zero() and poly_sign_definite(Q[a][b].d) != 1:
                info['status'] = f'DENOMINATOR NOT POSITIVE at {names[a]},{names[b]}'
                return False, info
    sig = RF(PL([0, 1]))
    kv = mlib.kernel_vectors(names, k, j, i, sig)
    kvR = [(lab, {nm: (cf if isinstance(cf, RF) else RF(flint.fmpq(cf))) for nm, cf in v.items()}) for lab, v in kv]
    if bscale == 1:
        for lab, v in kvR:
            for a in range(n):
                s = RF(0)
                for nm, cf in v.items():
                    s = s + Q[a][names.index(nm)] * cf
                if not s.is_zero():
                    info['status'] = f'KERNEL VECTOR {lab} NOT ANNIHILATED'
                    return False, info
    S, how = choose_gauge(names, kvR)
    if S is None:
        info['status'] = 'NO GAUGE'
        return False, info
    keep = [c for c in range(n) if c not in S]
    info['kernel'] = [lab for lab, _ in kv]
    info['gauge'] = [names[c] for c in S]
    info['gauge_det'] = how
    M = [[Q[a][b] for b in keep] for a in keep]
    info['form'] = [[(list(map(str, x.n.coeffs())), list(map(str, x.d.coeffs()))) for x in row] for row in M]
    pivs = []
    nk = len(keep)
    maxdeg = 0
    for t in range(nk):
        piv = M[t][t]
        if piv.is_zero():
            info['status'] = f'ZERO PIVOT at {names[keep[t]]}'
            return False, info
        ok, how2, wit = positive_on_halfline(piv)
        maxdeg = max(maxdeg, piv.n.degree(), piv.d.degree())
        if not ok:
            info['status'] = f'PIVOT NOT POSITIVE at {names[keep[t]]}: {how2}; witness sigma = {wit}'
            info['witness'] = wit
            info['failed_at'] = t
            return False, info
        pivs.append((names[keep[t]], how2))
        for a in range(t + 1, nk):
            if M[a][t].is_zero():
                continue
            f = M[a][t] / piv
            for b in range(t + 1, nk):
                if not M[t][b].is_zero():
                    M[a][b] = M[a][b] - f * M[t][b]
    info['status'] = 'PD on complement'
    info['pivots'] = pivs
    info['nkernel'] = len(kv)
    info['ngauged'] = nk
    info['maxdeg'] = maxdeg
    return True, info


if __name__ == '__main__':
    kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
    bscale = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1)
    save = sys.argv[4] if len(sys.argv) > 4 else None
    T0 = time.time()
    allok = True
    nforms = 0
    saved = {}
    nonneg_all = True
    failed_k = []
    for k in range(kmin, kmax + 1):
        t0 = time.time()
        cnt = 0
        md = 0
        kfail = None
        for j in range(-1, k + 1):
            for i in range(0, (k - 1) // 2 + 1):
                ok, info = certify(k, j, i, bscale=bscale)
                if info['status'] == 'empty':
                    continue
                cnt += 1
                nforms += 1
                if ok:
                    md = max(md, info['maxdeg'])
                    if any('Sturm' in h or 'negative coefficients' not in h for _, h in info['pivots']):
                        nonneg_all = False
                    if save:
                        saved[(k, j, i)] = info
                else:
                    allok = False
                    if bscale == 1:
                        print(f'  k={k} j={j} i={i}: {info["status"]}', flush=True)
                    if bscale != 1:
                        kfail = kfail or []
                        kfail.append((j, i, str(info.get('witness'))))
                        break
        if kfail is not None:
            failed_k.append(k)
            print(f'  k={k}: failing (j, i, sigma): {kfail}', flush=True)
        print(f'k={k}: {cnt} reduced forms, largest pivot degree {md}, {time.time() - t0:.1f}s', flush=True)
    if save:
        with open(save, 'wb') as fh:
            pickle.dump(saved, fh)
    if bscale == 1:
        print(f'{nforms} forms; every pivot with nonnegative coefficients: {nonneg_all}')
        print('PASS' if allok else 'FAIL', f'{time.time() - T0:.0f}s')
    else:
        nk = kmax - kmin + 1
        print(f'negative control bscale = {bscale}: some form fails for {len(failed_k)} of the {nk} values of k')
        print('FAIL (expected: the scaled block is not positive semidefinite)' if len(failed_k) == nk else
              'PASS (not every k failed)', f'{time.time() - T0:.0f}s')
