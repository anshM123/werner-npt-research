"""Fixed-k certificates for the type I blocks with a letter 1 and one marked copy.

For each k in [KMIN, KMAX] and every (j, tau, case, mode, i) of tlib.all_forms(k) (3 <= j <= k - 1, 1 <= tau <= j,
marked copy in T or in U, T'-mode sym or anti, every irrep index i with a nonzero twisted space):
  1. Q over Q(sigma) (tlib.reduced, symbolic mode); Q symmetric; every entry denominator > 0 on (0, inf);
  2. Q kappa = 0 identically for every known kernel vector kappa (tlib.kernel_vectors: the twisted parts of p (x) H_j);
  3. gauge: a set S of coordinates, one per kernel vector, with det(kappa_l[S]) a nonzero rational function of sigma,
     sign-definite on (0, inf); then the coordinate subspace W on the other coordinates is a complement of span(kappa)
     for every sigma > 0, and Q >= 0 with kernel exactly span(kappa) iff Q restricted to W is positive definite;
  4. symmetric elimination of Q on W in the order of tlib.NAMES; every pivot is > 0 on (0, inf)
     (signtest.positive_on_halfline: nonnegative coefficients, else a witness search, else Sturm).  Then every leading
     principal minor of Q on W is > 0 there (Sylvester).
Usage: python cert_t.py KMIN KMAX [bscale] [save.pkl] [JMIN JMAX]
With bscale != 1 (negative control) the kernel identities are not required; the run reports, per k, the forms whose
certificate fails, with a witness sigma, and FAIL if every k has a failing form (expected for 1001/1000)."""
import sys, time, itertools, pickle
from fractions import Fraction as Fr
import tlib
from ufield import RF
from signtest import positive_on_halfline
import flint

PL = flint.fmpq_poly


def poly_sign_definite(p):
    f = RF(p)
    ok, how, wit = positive_on_halfline(f)
    if ok:
        return 1
    ok, how, wit = positive_on_halfline(RF(-p))
    if ok:
        return -1
    return 0


def det_poly(M):
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


PREF = ['CM', 'C', 'O0', 'O1', 'C0', 'C1']


def choose_gauge(names, kvecs):
    r = len(kvecs)
    if r == 0:
        return [], 'none'
    rows = [[RF(flint.fmpq(v.get(nm, 0))) for nm in names] for _, v in kvecs]
    pref = sorted(range(len(names)), key=lambda c: (PREF.index(names[c]) if names[c] in PREF else len(PREF), c))
    for S in itertools.combinations(pref, r):
        sub = [[row[c] for c in S] for row in rows]
        d = det_poly(sub)
        if d.is_zero():
            continue
        if poly_sign_definite(d.n) != 0 and poly_sign_definite(d.d) != 0:
            return list(S), f'det {d}'
    return None, 'no gauge'


def certify(k, j, tau, case, mode, i, bscale=1):
    C = tlib.Calc(k, j, tau, i, case, mode, sigma=None, bscale=bscale)
    names, norms, Q = tlib.reduced(C)
    info = {'k': k, 'j': j, 'tau': tau, 'case': case, 'mode': mode, 'i': i, 'names': names}
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
    kv = tlib.kernel_vectors(C, names)
    kvR = [(lab, {nm: RF(flint.fmpq(cf)) for nm, cf in v.items()}) for lab, v in kv]
    if bscale == 1:
        for lab, v in kvR:
            for a in range(n):
                s = RF(0)
                for nm, cf in v.items():
                    s = s + Q[a][names.index(nm)] * cf
                if not s.is_zero():
                    info['status'] = f'KERNEL VECTOR {lab} NOT ANNIHILATED'
                    return False, info
    S, how = choose_gauge(names, kv)
    if S is None:
        info['status'] = 'NO GAUGE'
        return False, info
    keep = [c for c in range(n) if c not in S]
    info['kernel'] = [lab for lab, _ in kv]
    info['kvecs'] = kv
    info['gauge'] = [names[c] for c in S]
    info['gauge_det'] = how
    info['keep'] = [names[c] for c in keep]
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
        # symmetric elimination: the Schur complement stays symmetric, so only b >= a is computed and mirrored
        for a in range(t + 1, nk):
            if M[t][a].is_zero():
                continue
            f = M[t][a] / piv
            for b in range(a, nk):
                if not M[t][b].is_zero():
                    M[a][b] = M[a][b] - f * M[t][b]
                    if b != a:
                        M[b][a] = M[a][b]
    info['status'] = 'PD on complement'
    info['pivots'] = pivs
    info['nkernel'] = len(kv)
    info['ngauged'] = nk
    info['maxdeg'] = maxdeg
    return True, info


if __name__ == '__main__':
    kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
    bscale = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1)
    save = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] != '-' else None
    jlo = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    jhi = int(sys.argv[6]) if len(sys.argv) > 6 else None
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
        kfail = []
        forms = tlib.all_forms(k, jlo, (k - 1) if jhi is None else min(jhi, k - 1))
        for (j, tau, case, mode, i) in forms:
            ok, info = certify(k, j, tau, case, mode, i, bscale=bscale)
            if info['status'] == 'empty':
                continue
            cnt += 1
            nforms += 1
            if ok:
                md = max(md, info['maxdeg'])
                if any('nonnegative coefficients' not in h for _, h in info['pivots']):
                    nonneg_all = False
                if save:
                    saved[(k, j, tau, case, mode, i)] = info
            else:
                allok = False
                if bscale == 1:
                    print(f'  k={k} j={j} tau={tau} {case} {mode} i={i}: {info["status"]}', flush=True)
                else:
                    kfail.append((j, tau, case, mode, i, str(info.get('witness'))))
        if bscale != 1:
            if kfail:
                failed_k.append(k)
            print(f'  k={k}: {len(kfail)} of {cnt} forms fail; first: {kfail[:3]}', flush=True)
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
