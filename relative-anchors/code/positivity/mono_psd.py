"""Test 1 (exact): does den * M_g have a representation  sum_e c_e(s) x_e x_e^T  with fixed rational vectors x_e and
weights c_e that are polynomials with nonnegative coefficients?  With den = lcm of the entry denominators and no extra
multiplier this holds iff every monomial-coefficient matrix C_m of den * M_g is positive semidefinite (C_m = sum of
the x_e x_e^T carrying that monomial; conversely an LDL^T of each C_m gives the x_e).  The same test is run with the
multipliers |s|^N (N = 1, 2) and, for small blocks, with every product of distinct den factors.
Also prints the sign profile of every entry (numerator and denominator), the factorisation of den, and the off-diagonal
sign pattern.  Usage: python mono_psd.py k [maxsize]"""
import sys, time, itertools
from fractions import Fraction as Fr
import flint
from psdtools import load_blocks, poly_matrix, coeff_mats, psd_exact, sign_profile

k = int(sys.argv[1])
maxsize = int(sys.argv[2]) if len(sys.argv) > 2 else 99
t0 = time.time()
blocks = load_blocks(k)
summary = []
vertex_counts = []
for blk in blocks:
    n = len(blk['Mg'])
    if n > maxsize:
        continue
    ctx = blk['ctx']
    Mg = blk['Mg']
    den, Nt = poly_matrix(blk)
    cst, fac = den.factor()
    facs_pos = all(all(c > 0 for c in f.coeffs()) for f, e in fac) and cst > 0
    diag = [sign_profile(Nt[a][a]) for a in range(n)]
    offd = {}
    for a in range(n):
        for c in range(a + 1, n):
            offd.setdefault(sign_profile(Nt[a][c]), []).append((a, c))
    print(f'== {blk["name"]} ({blk["kind"]}), gauge-fixed size {n}', flush=True)
    print(f'   den = {cst} * ' + ' * '.join(f'[deg {f.total_degree()}, {len(f)} t, {sign_profile(f)}]^{e}' for f, e in fac)
          + f'; den factors positive: {facs_pos}')
    print(f'   diagonal of den*M sign profiles: {"".join(d if d != "mixed" else "?" for d in diag)} '
          f'(? = mixed);  off-diagonal: ' + ', '.join(f'{key}: {len(v)}' for key, v in sorted(offd.items())))
    res = []
    for tag, mult in [('q=1', None), ('q=|s|', None), ('q=|s|^2', None)]:
        S1 = sum((g for g in ctx.gens()), ctx.from_dict({}))
        if tag == 'q=|s|':
            mult = S1
        elif tag == 'q=|s|^2':
            mult = S1 * S1
        if k == 4 and tag != 'q=1':
            # dehomogenised (s0 = 1): |s| = 1 + s1 + s2 + s3
            one = ctx.from_dict({(0,) * len(blk['names']): 1})
            S1 = one + S1
            mult = S1 if tag == 'q=|s|' else S1 * S1
        _, Ntm = poly_matrix(blk, mult)
        C = coeff_mats(Ntm)
        nfail = 0
        first = None
        for m in sorted(C):
            okm, wit = psd_exact(C[m])
            if not okm:
                nfail += 1
                if first is None:
                    first = (m, wit)
        res.append((tag, len(C), nfail))
        nvert = 0
        if tag == 'q=1' and nfail:
            # vertex obstruction: a non-PSD coefficient matrix at a vertex of the Newton polytope of den*M_g excludes
            # every fixed-vector representation with subtraction-free weights (exact certificate: direction l)
            supp = list(C)
            nvv = len(supp[0])
            dirs = [l for l in itertools.product(range(-4, 5), repeat=nvv) if any(l)]
            for m in supp:
                if psd_exact(C[m])[0]:
                    continue
                for l in dirs:
                    vm = sum(x * y for x, y in zip(l, m))
                    if all(sum(x * y for x, y in zip(l, m2)) < vm for m2 in supp if m2 != m):
                        nvert += 1
                        break
            vertex_counts.append((blk['name'], nvert))
        msg = f'   {tag}: {len(C)} monomials, {nfail} with C_m not PSD' + (f' ({nvert} of them at vertices of the Newton polytope)' if tag == 'q=1' and nfail else '')
        if first is not None:
            m, wit = first
            sup = [(blk['labels'][i], str(v)) for i, v in enumerate(wit) if v != 0]
            msg += f'; e.g. monomial {m}: witness {sup[:6]}{"..." if len(sup) > 6 else ""}'
        print(msg, flush=True)
    summary.append((blk['name'], blk['kind'], n, facs_pos, diag.count('mixed'), len(offd.get('mixed', [])), res))
    print(f'   [{time.time() - t0:.1f}s]', flush=True)
print()
print('SUMMARY (block, kind, size, den factors positive, #diag mixed, #offdiag mixed, [(multiplier, #monomials, #non-PSD)])')
for row in summary:
    print('  ', row)
npass = sum(1 for row in summary if row[6][0][2] == 0)
print(f'VERTEX k = {k}: blocks with a non-PSD coefficient matrix at a vertex (no fixed-vector representation with '
      f'subtraction-free weights): {[nm for nm, c in vertex_counts if c > 0]}')
print(f'RESULT k = {k}: {npass} of {len(summary)} blocks have every coefficient matrix of den*M_g PSD (q = 1); '
      f'passing blocks: {[row[0] for row in summary if row[6][0][2] == 0]}')
print(f'time {time.time() - t0:.1f}s')
print('PASS')
