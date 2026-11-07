"""Coverage bookkeeping of the type I forms with a letter 1 at one marked copy, for every k <= KCHK.

Every (j, tau, case, mode, i) with 3 <= j <= k - 1, 1 <= tau <= j, the marked copy in T (MT) or in U (MU), T'-mode
sym or anti and i = 0..floor(u'/2) is enumerated through its sizes (q0, t, i, tau'), with k = tau' + q0 + t + 2i + 1.
Its twisted space is nonzero iff some twisted vector is present; presence is computed by a closed form of
tlib.norm_of, which is compared with tlib.norm_of itself for every form with k <= KCMP.  A form with j = k - 1 is
outside the scope here (Corollary M1) and only counted.  A form with 3 <= j <= k - 2 is covered
  - by Theorem I1 if k <= KFIX (29: the fixed-k certificates of ../typeI);
  - by a family if some family certified 'PD' in the pickles contains it (each size equal to the fixed value, or at
    least the minimum of a symbolic size).
Pickles of q3cert.py (key 'results3') certify only the three-row part Q3 of a family: they cover no form, and for the
forms left open the script counts those whose canonical family has Q3 certified (or no split block).
Checked: every nonempty form has q0 >= -1, t >= -1 (MT) or t >= -2 (MU), i >= 0 and tau' >= 0 (MT) or >= 1 (MU), so
it lies in exactly one family of famlist.py (its canonical family) and that family is listed there when k >= 30; no
family reported empty contains a nonempty form.  Printed: per (case, mode) the forms with 30 <= k <= KCHK covered by
families and not covered, the latter by the number of symbolic sizes of their canonical family; the smallest k of an
uncovered form; and, for every k in a list, the counts.
python coverage_st3.py KCHK KFIX KCMP PKL [PKL ...]"""
import sys, pickle
from collections import Counter, defaultdict
import tlib
import tfamlib as T
import famlist


def cb(n, r):
    return 0 <= r <= n


def present(case, mode, q0, t, i, tp):
    m = q0 + t
    gok = (t + i >= 1) if case == 'MT' else (t + i >= 0)   # j + 1 <= k
    if tp < 0 or m < 0 or i < 0:
        return []
    if mode == 'anti' and tp < 2:
        return []
    out = []
    if case == 'MT':
        if mode == 'sym':
            for x in (0, 1):
                if m >= 1 and cb(m - 1, q0 - x):
                    out.append('A%d' % x)
            for x in (0, 1):
                for y in (0, 1):
                    if i >= 1 and cb(m, q0 + 1 - x - y):
                        out.append('B%d%d' % (x, y))
            if tp >= 1 and cb(m, q0):
                out.append('C')
            for z in (0, 1):
                if tp >= 1 and cb(m, q0 + 1 - z):
                    out.append('D%d' % z)
            if cb(m, q0):
                out.append('CM')
            for z in (0, 1):
                if cb(m, q0 + 1 - z):
                    out.append('DM%d' % z)
            if cb(m, q0 + 1) and gok:
                out.append('G')
        else:
            if cb(m, q0):
                out.append('Ca')
            for z in (0, 1):
                if cb(m, q0 + 1 - z):
                    out.append('Da%d' % z)
        return out
    if mode == 'sym':
        for y0 in (0, 1):
            if cb(m, q0 - y0):
                out.append('O%d' % y0)
        for x in (0, 1):
            for y0 in (0, 1):
                if m >= 1 and cb(m - 1, q0 - x - y0):
                    out.append('A%d%d' % (x, y0))
        for x in (0, 1):
            for y in (0, 1):
                for y0 in (0, 1):
                    if i >= 1 and cb(m, q0 + 1 - x - y - y0):
                        out.append('B%d%d%d' % (x, y, y0))
        for y0 in (0, 1):
            if tp >= 1 and cb(m, q0 - y0):
                out.append('C%d' % y0)
        for z in (0, 1):
            for y0 in (0, 1):
                if tp >= 1 and cb(m, q0 + 1 - z - y0):
                    out.append('D%d%d' % (z, y0))
        for y0 in (0, 1):
            if cb(m, q0 + 1 - y0) and gok:
                out.append('G%d' % y0)
        return out
    for y0 in (0, 1):
        if cb(m, q0 - y0):
            out.append('Ca%d' % y0)
    for z in (0, 1):
        for y0 in (0, 1):
            if cb(m, q0 + 1 - z - y0):
                out.append('Da%d%d' % (z, y0))
    return out


def forms(k):
    """every (j, tau, case, mode, i, sizes) with 3 <= j <= k - 1, 1 <= tau <= j"""
    for j in range(3, k):
        for tau in range(1, j + 1):
            for case in ('MT', 'MU'):
                up = k - tau if case == 'MT' else k - tau - 1
                if up < 0:
                    continue
                tp = tau - 1 if case == 'MT' else tau
                for mode in ('sym', 'anti'):
                    for i in range(0, up // 2 + 1):
                        m = up - 2 * i
                        q0 = j - tau - i
                        t = m - q0
                        yield j, tau, case, mode, i, (q0, t, i, tp)


def main():
    kchk, kfix, kcmp = (int(x) for x in sys.argv[1:4])
    fams = {}
    empty = []
    q3 = {}
    for p in sys.argv[4:]:
        d = pickle.load(open(p, 'rb'))
        if 'results3' in d:
            # certificates of the three-row parts only (q3cert.py): not counted as covering a form
            for key, r in d['results3'].items():
                if r.startswith('Q3 positive') or r.startswith('no split'):
                    q3[key[:2] + tuple(None if isinstance(x, tuple) else x for x in key[2:])] = r
            continue
        for key, r in d['results'].items():
            if r == 'empty':
                empty.append(key)
            elif r.startswith('PD'):
                fams[key] = r
    proved = defaultdict(list)
    for key in fams:
        proved[key[:2]].append(key[2:])
    emp = defaultdict(list)
    for key in empty:
        emp[key[:2]].append(key[2:])

    def inside(spec, sizes):
        for s, v in zip(spec, sizes):
            if isinstance(s, tuple):
                if v < s[1]:
                    return False
            elif s != v:
                return False
        return True

    ok = True
    # 1. the closed form of presence against tlib.norm_of
    ncmp = 0
    for k in range(4, kcmp + 1):
        for j, tau, case, mode, i, sz in forms(k):
            C = tlib.Calc(k, j, tau, i, case, mode, sigma=1)
            a = [nm for nm in tlib.NAMES[(case, mode)] if C.valid() and tlib.norm_of(C, nm) != 0]
            b = present(case, mode, *sz)
            ncmp += 1
            if a != b:
                ok = False
                print('PRESENCE MISMATCH', k, j, tau, case, mode, i, a, b)
    print(f'closed form of the presence of twisted vectors equals tlib.norm_of on all {ncmp} forms with k <= {kcmp}',
          flush=True)
    listed = {(c, mo): set(famlist.families(c, mo)) for c in ('MT', 'MU') for mo in ('sym', 'anti')}
    other_mins = {cm: [s for s in specs if any(isinstance(x, tuple) and x[1] != T.SMIN[nm]
                                               for x, nm in zip(s, T.SIZES))] for cm, specs in proved.items()}
    default = {cm: {tuple(None if isinstance(x, tuple) else x for x in s) for s in specs
                    if s not in other_mins[cm]} for cm, specs in proved.items()}
    CAP = 8
    pcache, fcache = {}, {}
    cnt = Counter()
    unc = Counter()
    q3cov = Counter()
    kmin_unc = {}
    per_k = {}
    for k in range(4, kchk + 1):
        row = Counter()
        for j, tau, case, mode, i, sz in forms(k):
            q0, t, _, tp = sz
            ck = (case, mode) + tuple(min(v, CAP) for v in sz)
            pres = pcache.get(ck)
            if pres is None:
                pres = pcache[ck] = bool(present(case, mode, *ck[2:]))
            if not pres:
                continue
            if q0 < -1 or t < (-1 if case == 'MT' else -2) or tp < (0 if case == 'MT' else 1):
                ok = False
                print('NONEMPTY FORM OUTSIDE THE SIZE RANGE', k, j, tau, case, mode, i, sz)
            canon = tuple(None if v >= T.SMIN[nm] else v for v, nm in zip(sz, T.SIZES))
            fk = (case, mode, canon)
            st = fcache.get(fk)
            if st is None:
                st = fcache[fk] = {'listed': canon in listed[(case, mode)],
                                   'proved': canon in default.get((case, mode), set()),
                                   'empty': [s for s in emp[(case, mode)]]}
            for spec in st['empty']:
                if inside(spec, sz):
                    ok = False
                    print('FAMILY REPORTED EMPTY CONTAINS A NONEMPTY FORM', spec, k, j, tau, case, mode, i)
            if j == k - 1:
                row['j=k-1 (Corollary M1)'] += 1
                continue
            if k >= 30 and not st['listed']:
                ok = False
                print('CANONICAL FAMILY NOT LISTED', k, j, tau, case, mode, i, sz, canon)
            if k <= kfix:
                row['Theorem I1'] += 1
                continue
            if st['proved'] or any(inside(spec, sz) for spec in other_mins.get((case, mode), [])):
                row['family'] += 1
                cnt[(case, mode, 'family')] += 1
                continue
            row['open'] += 1
            ns = sum(v is None for v in canon)
            unc[(case, mode, ns)] += 1
            if mode == 'anti' or (case, mode) + canon in q3:
                q3cov[(case, mode)] += 1
            kmin_unc[(case, mode)] = min(kmin_unc.get((case, mode), 10 ** 9), k)
        per_k[k] = row
    for k in (8, 29, 30, 40, 60, 80, 100, 120):
        if k in per_k:
            print(f'k={k}: ' + ', '.join(f'{a} {b}' for a, b in sorted(per_k[k].items())), flush=True)
    tot = Counter()
    for r in per_k.values():
        tot.update(r)
    print(f'k = 4..{kchk}: ' + ', '.join(f'{a} {b}' for a, b in sorted(tot.items())))
    for c in ('MT', 'MU'):
        for mo in ('sym', 'anti'):
            u = {ns: unc[(c, mo, ns)] for ns in range(1, 5) if unc[(c, mo, ns)]}
            print(f'{c} {mo}: 30 <= k <= {kchk}: {cnt[(c, mo, "family")]} forms in proved families; not covered: '
                  f'{sum(u.values())} (by symbolic sizes of the canonical family: {u}); smallest k not covered: '
                  f'{kmin_unc.get((c, mo), "none")}')
    if q3:
        for c in ('MT', 'MU'):
            nu_ = sum(unc[(c, 'sym', ns)] for ns in range(1, 5))
            print(f'{c} sym: of the {nu_} forms not covered, {q3cov[(c, "sym")]} lie in a family whose three-row part '
                  f'Q3 is certified (or has no split block)')
    print(f'families read: {len(fams)} proved, {len(empty)} empty, {len(q3)} with the three-row part certified')
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
