"""Run the certificate of cert.py for all 16 families at several regimes, in one process.

python batch.py REGIME [REGIME ...]    a regime is U0,U1 with U = an integer or 'o+' (e.g. 9+,2  6+,6+)
Writes one line per (regime, family) with PASS/FAIL and the sizes, the forms to logs/forms_<U0>_<U1>.pkl, and PASS at
the end iff every certificate passed."""
import sys, time, pickle
import r2v
from cert import build, signcheck, det

t0 = time.time()


def certify(fam, s0, s1):
    K, ctx, F, bases, red = build(fam, s0, s1, log=lambda m: None)
    one = K(1)
    nK, nW, nQ = (len(b) for b in bases)
    info = dict(nK=nK, nW=nW, nQ=nQ)
    if nW == 0:
        return True, info, None
    ok = True
    for gname in ('GP', 'GQ'):
        if gname in red:
            sg, _ = signcheck(det(red[gname], one))
            ok &= sg == 1
    Sf = red['Sfull']
    n = len(Sf)
    ok &= all(Sf[i][j] == Sf[j][i] for i in range(n) for j in range(i + 1, n))
    ok &= all(Sf[i][j] == 0 for i in range(nK) for j in range(n))
    if not ok:
        info['why'] = 'Gram/symmetry/kernel'
        return False, info, None
    A = [list(r) for r in red['S']]
    m = len(A)
    big = (0, 0)
    for c in range(m):
        p = A[c][c]
        sg, desc = signcheck(p)
        big = max(big, (len(p.n), p.n.total_degree()))
        if sg != 1:
            info['why'] = f'pivot {c + 1}: {desc}'
            return False, info, None
        for i in range(c + 1, m):
            if not (A[i][c] == 0):
                f = A[i][c] / p
                for j in range(c + 1, m):
                    A[i][j] = A[i][j] - f * A[c][j]
    info['largest pivot numerator (terms, degree)'] = big
    form = dict(names=K.ctx.names(), basis=red['names'], S=[[(str(x.n), str(x.d)) for x in row] for row in red['S']])
    return True, info, form


if __name__ == '__main__':
    allok = True
    for reg in sys.argv[1:]:
        s0, s1 = reg.split(',')
        forms = {}
        for F in r2v.families():
            ok, info, form = certify(F.name, s0, s1)
            allok &= ok
            forms[F.name] = form
            print(f'({s0}, {s1}) {F.name}: {"PASS" if ok else "FAIL"} {info} [{time.time() - t0:.1f}s]', flush=True)
        with open(f'logs/forms_{s0}_{s1}.pkl', 'wb') as fh:
            pickle.dump(forms, fh)
    print(f'time {time.time() - t0:.1f}s')
    print('PASS' if allok else 'FAIL')
