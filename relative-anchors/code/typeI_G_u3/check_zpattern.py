"""Exact check of the sign pattern of G_z for m = 1 (u = 1, 2, 3; closed forms closed_u{u}_m1.pkl).

Claim checked: G_z has one common denominator with only positive coefficients, every diagonal numerator has only
positive coefficients and every off-diagonal numerator has only negative coefficients.  So for m = 1 and u <= 3, G_z is
a symmetric Z-matrix on the open orthant, and G > 0 is equivalent to G_z being a nonsingular M-matrix (some v > 0 with
G_z v > 0).  This is an observation that may point to a proof for every u when m = 1 (NOTES.md, Section 4b); the
theorem for u = 3 does not use it.  For m = 2 the pattern fails (entries of G_tt' with mixed coefficients); printed."""
import pickle, flint, os


def sg(p):
    c = list(p.to_dict().values())
    return '+' if c and all(x > 0 for x in c) else ('-' if c and all(x < 0 for x in c) else ('0' if not c else '?'))


ok = True
for (u, m) in [(1, 1), (2, 1), (3, 1), (1, 2), (2, 2)]:
    fn = f'closed_u{u}_m{m}.pkl'
    D = pickle.load(open(fn, 'rb'))
    ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
    G = D['G']
    n = len(G)
    nums = [[ctx.from_dict(G[i][j][0]) for j in range(n)] for i in range(n)]
    dens = [[ctx.from_dict(G[i][j][1]) for j in range(n)] for i in range(n)]
    pat = [[sg(nums[i][j]) + sg(dens[i][j]) for j in range(n)] for i in range(n)]
    z = all(pat[i][j] == ('++' if i == j else '-+') for i in range(n) for j in range(n))
    print(f'(u, m) = ({u}, {m}): Z-matrix pattern (diagonal ++, off-diagonal -+): {z}')
    for row in pat:
        print('    ' + ' '.join(row))
    if m == 1:
        ok &= z
print('PASS' if ok else 'FAIL')
