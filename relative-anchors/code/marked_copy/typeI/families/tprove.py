"""Family certificates for the type I blocks with a letter 1 at one marked copy (tfamlib.py).

A family fixes each of q0, t, i, tau' or leaves it symbolic with a minimum (tfamlib.SMIN unless changed); symbolic
sizes are shifted to their minimum (q0 = qmin + x, ...), sigma > 0 is symbolic.  For a family:
  1. the form Q over Q(sizes, sigma) (tfamlib.reduced); Q symmetric;
  2. sanity: N0 times Q equals the exact numeric form of tlib.py (which asserts its Newton window) at every sample
     point (each symbolic size at min, min + 1, min + 3) and sigma in {2/7, 5}, with the same twisted vectors;
  3. every intermediate denominator of the calculus (tfamlib C.dens) and every entry denominator is one-signed after
     the shifts (all coefficients of one sign, a nonzero part free of the size variables), so it does not vanish at
     any member and any sigma > 0;
  4. the known kernel vectors (tlib.kernel_vectors) are annihilated identically;
  5. gauge: one deleted coordinate per kernel vector, the determinant of the kernel vectors there a nonzero constant;
  6. greedy symmetric elimination on the other coordinates: every pivot a ratio of one-signed polynomials (after the
     shifts) of the same sign, positive.
Then at every member and every sigma > 0 the reduced form is >= 0 with kernel exactly the span of the known kernel
vectors (Sylvester on the gauge-fixed form, whose leading principal minors in the order taken are products of pivots).
Usage: python tprove.py OUT.pkl FAM [FAM ...] [--bscale=B] [--tlimit=S] [--smin=q0:4,t:4,i:2,tp:2]
FAM = case/mode/q0,t,i,tp with 'n' for a symbolic size, e.g. MT/sym/n,2,1,n"""
import sys, time, itertools, pickle
from fractions import Fraction as Fr
import flint
import tfamlib as T
import tlib

PREF = ['CM', 'O0', 'O1', 'C', 'C0', 'C1']


def shifted(P, C):
    gens = C.K.ctx.gens()
    args = [g + C.smin[nm] if nm in T.SIZES else g for g, nm in zip(gens, C.syms)]
    return P.compose(*args)


def one_signed(P, C):
    """+1 / -1 if after the shifts all coefficients have that sign and the part free of the size variables is not
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


def fam_label(case, mode, fam, smin=None):
    smin = smin or T.SMIN
    return f'{case} {mode} [' + ', '.join(f'{nm}>={smin[nm]}' if v is None else f'{nm}={v}'
                                       for v, nm in zip(fam, T.SIZES)) + ']'


def sample_points(C, extra=()):
    sz = [nm for nm in C.syms if nm != 's']
    opts = [[C.smin[nm] + o for o in ((0, 1, 3) + tuple(extra))] for nm in sz]
    for sample in itertools.product(*opts):
        vals = {nm: C.fixed[nm] for nm in T.SIZES}
        vals.update(dict(zip(sz, sample)))
        yield [vals[nm] for nm in sz], vals


def sanity(C, names, Q, bscale, sigmas=(Fr(2, 7), Fr(5)), points=None):
    """N0 Q equals tlib's numeric form at the sample points; returns (None or an error, number of points)"""
    npts = 0
    for pt, vals in (points if points is not None else sample_points(C)):
        k, j, tau = C.member(vals)
        N0 = T.norm0(C, vals)
        for sg in sigmas:
            Cn = tlib.Calc(k, j, tau, vals['i'], C.case, C.mode, sigma=sg, bscale=bscale)
            n2, _, Q2 = tlib.reduced(Cn)
            if n2 != names:
                return f'SANITY-NAMES {names} {n2} at {(k, j, tau, vals["i"])}', npts
            ptv = list(pt) + [sg]
            for a in range(len(names)):
                for b in range(len(names)):
                    if evalF(Q[a][b], ptv) * N0 != Fr(int(Q2[a][b].p), int(Q2[a][b].q)):
                        return f'SANITY-VALUES at {(k, j, tau, vals["i"], sg)} {names[a]},{names[b]}', npts
            npts += 1
    return None, npts


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


def tsz(x):
    return 0 if x == 0 else len(x.n) + len(x.d)


def gauge_for(C, names, kv):
    r = len(kv)
    if r == 0:
        return []
    pref = sorted(range(len(names)), key=lambda c: (PREF.index(names[c]) if names[c] in PREF else len(PREF), c))
    for S in itertools.combinations(pref, r):
        sub = [[v.get(names[c], C.zero) for c in S] for _, v in kv]
        d = det_small(sub, C)
        if d != 0 and d.n.is_constant() and d.d.is_constant():
            return list(S)
    return None


def eliminate(M, act, C, t0, tlimit=None, verbose=False, names=None):
    act = list(act)
    info = {'order': [], 'maxterms': 0, 'pivots': []}
    while act:
        if tlimit is not None and time.time() - t0 > tlimit:
            return 'TIME', info
        cands = sorted(act, key=lambda z: tsz(M[z, z]) if M[z, z] != 0 else 10 ** 12)
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
            nm = names[x] if names else x
            return (f'PIVOT-SIGN-UNPROVED (smallest {nm}: num {len(M[x, x].n)} terms, den {len(M[x, x].d)})'
                    if M[x, x] != 0 else f'PIVOT-SIGN-UNPROVED (zero diagonal at {nm})'), info
        x = chosen
        piv = M[x, x]
        info['maxterms'] = max(info['maxterms'], len(piv.n), len(piv.d))
        info['order'].append(names[x] if names else x)
        info['pivots'].append((names[x] if names else x, len(piv.n), len(piv.d), round(time.time() - t0, 1)))
        if verbose:
            print(f'    pivot {names[x] if names else x}: num {len(piv.n)} terms, den {len(piv.d)} terms '
                  f'({time.time()-t0:.1f}s)', flush=True)
        act.remove(x)
        for ia, a in enumerate(act):
            if M[a, x] == 0:
                continue
            f = M[a, x] / piv
            for b in act[ia:]:
                if M[x, b] != 0:
                    M[a, b] = M[a, b] - f * M[x, b]
                    M[b, a] = M[a, b]
    return 'OK', info


def prove(case, mode, fam, bscale=Fr(1), tlimit=None, verbose=False, smin=None, do_sanity=True):
    """the certificate of the module docstring.  In the anti mode a family with i and tau' symbolic is reduced to
    i = 2 by the structure lemma (structure.py: with tau' = nu - 2i its form is free of i, so the member (i, tau')
    has the form of the member (2, tau' + 2i - 4), and tau' + 2i - 4 >= tau' >= 2): the lemma is checked exactly on
    the form with i and tau' symbolic, then the family with i = 2 is certified."""
    if mode == 'anti' and fam[2] is None and fam[3] is None and (smin or T.SMIN)['i'] == 2:
        import structure
        t0 = time.time()
        try:
            good, worst, dmax, nv = structure.check(case, mode, fam)
        except T.NonUniform as ex:
            return f'NONUNIFORM {ex}', time.time() - t0, {}
        if not good:
            return 'STRUCTURE LEMMA FAILS', time.time() - t0, {}
        res, dt, info = prove(case, mode, (fam[0], fam[1], 2, None), bscale, tlimit, verbose, smin, do_sanity)
        info['structure'] = f'i-free with tau1 = nu - 2i ({nv} vectors); certified at i = 2'
        return res, time.time() - t0, info
    t0 = time.time()
    q0, t, i, tp = fam
    try:
        C = T.FamCalc(case, mode, q0, t, i, tp, bscale=bscale, smin=smin)
        if not C.valid():
            return 'empty', 0, {}
        names, Q = T.reduced(C)
    except T.NonUniform as ex:
        return f'NONUNIFORM {ex}', time.time() - t0, {}
    if not names:
        return 'empty', time.time() - t0, {}
    n = len(names)
    big = max(tsz(Q[a][b]) for a in range(n) for b in range(n))
    info = {'names': names, 'syms': C.syms, 'stats': {'sectors': len(C.stats['sectors']),
                                                     'd': (C.stats['dmin'], C.stats['dmax'])},
            'form_terms': big, 'form_time': round(time.time() - t0, 1), 'ndens': len(C.dens)}
    if verbose:
        print(f'  form computed ({time.time()-t0:.1f}s), {n} vectors, largest entry {big} terms, '
              f'{len(C.dens)} intermediate denominators', flush=True)
    for a in range(n):
        for b in range(a + 1, n):
            if Q[a][b] != Q[b][a]:
                return f'NOT SYMMETRIC at {names[a]},{names[b]}', time.time() - t0, info
    if do_sanity:
        bad, npts = sanity(C, names, Q, bscale)
        if bad:
            return bad, time.time() - t0, info
        info['sanity_points'] = npts
    for key, (x, why) in C.dens.items():
        if one_signed(x.n, C) == 0 or not x.d.is_constant():
            return f'INTERMEDIATE DENOMINATOR NOT ONE-SIGNED ({why})', time.time() - t0, info
    for a in range(n):
        for b in range(n):
            if Q[a][b] != 0 and one_signed(Q[a][b].d, C) == 0:
                return f'ENTRY DENOMINATOR NOT ONE-SIGNED at {names[a]},{names[b]}', time.time() - t0, info
    kv = T.kernel_vectors(C, names)
    if bscale != 1:
        info['kernel_annihilated'] = all(
            sum((Q[a][names.index(nm)] * cf for nm, cf in v.items()), C.zero) == 0 for lab, v in kv for a in range(n))
        # <kappa, Q kappa> at the minimal member and sigma = 1 (negative: the scaled form is not >= 0 there)
        mp = C.minpoint()
        ptv = [mp[nm] for nm in C.syms if nm != 's'] + [1]
        qk = []
        for lab, v in kv:
            s_ = C.zero
            for na, ca in v.items():
                for nb, cb_ in v.items():
                    s_ = s_ + ca * cb_ * Q[names.index(na)][names.index(nb)]
            qk.append(evalF(s_, ptv))
        info['kQk'] = qk
    if bscale == 1:
        for lab, v in kv:
            for a in range(n):
                s = C.zero
                for nm, cf in v.items():
                    s = s + Q[a][names.index(nm)] * cf
                if s != 0:
                    return f'KERNEL VECTOR {lab} NOT ANNIHILATED', time.time() - t0, info
    S = gauge_for(C, names, kv)
    if S is None:
        return 'NO GAUGE', time.time() - t0, info
    info['kernel'] = [lab for lab, _ in kv]
    info['gauge'] = [names[c] for c in S]
    keep = [c for c in range(n) if c not in S]
    M = {(a, b): Q[a][b] for a in keep for b in keep}
    res, einfo = eliminate(M, keep, C, t0, tlimit, verbose, names)
    info.update(einfo)
    if res != 'OK':
        return res, time.time() - t0, info
    return f'PD on complement (pivots {len(keep)}, kernel {len(kv)})', time.time() - t0, info


def parse_fam(s):
    case, mode, sizes = s.split('/')
    vals = tuple(None if x == 'n' else int(x) for x in sizes.split(','))
    return case, mode, vals


if __name__ == '__main__':
    opts = {'bscale': '1', 'tlimit': '', 'smin': '', 'allow-time': '0'}
    args = []
    for a in sys.argv[1:]:
        if a.startswith('--'):
            k_, v_ = a[2:].split('=')
            opts[k_] = v_
        else:
            args.append(a)
    out = args[0]
    fams = []
    for a in args[1:]:
        if a.startswith('@'):
            fams += [parse_fam(s) for s in open(a[1:]).read().split()]
        else:
            fams.append(parse_fam(a))
    bscale = Fr(opts['bscale'])
    tlimit = float(opts['tlimit']) if opts['tlimit'] else None
    smin = dict(T.SMIN)
    if opts['smin']:
        for it in opts['smin'].split(','):
            nm, v = it.split(':')
            smin[nm] = int(v)
    T0 = time.time()
    results, infos = {}, {}
    allok = True
    npr = nem = 0
    for case, mode, fam in fams:
        res, dt, info = prove(case, mode, fam, bscale=bscale, tlimit=tlimit, smin=smin)
        key = (case, mode) + tuple(('>=', smin[nm]) if v is None else v for v, nm in zip(fam, T.SIZES))
        results[key] = res
        infos[key] = info
        ok = res.startswith('PD') or res == 'empty' or (res == 'TIME' and opts['allow-time'] == '1')
        allok &= ok
        npr += res.startswith('PD')
        nem += res == 'empty'
        extra = ''
        if info.get('names'):
            extra = (f'  [{len(info["names"])} vectors, gauge {info.get("gauge")}, largest entry {info.get("form_terms")}'
                     f' terms, largest pivot {info.get("maxterms")} terms]')
        if bscale != 1 and 'kernel_annihilated' in info:
            extra += (f'  [kernel vectors annihilated with B scaled: {info["kernel_annihilated"]}; <kappa, Q kappa> '
                      f'at the minimal member, sigma = 1: {[float(x) for x in info.get("kQk", [])]}]')
        print(f'{fam_label(case, mode, fam, smin)} {res}  ({dt:.1f}s){extra}', flush=True)
    if bscale == 1:
        with open(out, 'wb') as fh:
            pickle.dump({'results': results, 'info': infos, 'smin': smin}, fh)
    nt = sum(1 for r in results.values() if r == 'TIME')
    print(f'{npr} families proved, {nem} empty, {nt} stopped by the time limit (open), of {len(fams)}; '
          f'{time.time() - T0:.0f}s')
    if bscale != 1:
        nf = len(fams) - npr - nem
        print(f'negative control bscale = {bscale}: the certificate fails for {nf} of the {len(fams) - nem} families')
        print('FAIL (expected: some certificate fails)' if nf > 0 else 'PASS (no certificate failed)')
    else:
        print('PASS' if allok else 'FAIL')
