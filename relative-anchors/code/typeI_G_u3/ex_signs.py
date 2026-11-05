"""Exploration: sign patterns of the entries of G_z (closed forms closed_u{u}_m{m}.pkl).  For each entry: are all
coefficients of the numerator of one sign?  Also the 'row excess' G_ii - sum_{j != i} |G_ij| when the off-diagonal
entries have one sign."""
import pickle, sys, flint
for (u, m) in [(1, 1), (2, 1), (1, 2), (2, 2), (3, 1)]:
    D = pickle.load(open(f'closed_u{u}_m{m}.pkl', 'rb'))
    ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
    G = D['G']
    n = len(G)
    pat = []
    for i in range(n):
        row = []
        for j in range(n):
            num = ctx.from_dict(G[i][j][0]); den = ctx.from_dict(G[i][j][1])
            cn = [c for c in num.to_dict().values()]; cd = [c for c in den.to_dict().values()]
            sn = '+' if all(c > 0 for c in cn) else ('-' if all(c < 0 for c in cn) else ('0' if not cn else '?'))
            sd = '+' if all(c > 0 for c in cd) else ('-' if all(c < 0 for c in cd) else '?')
            row.append(sn + sd)
        pat.append(row)
    print(f'(u, m) = ({u}, {m}): sign pattern (numerator sign, denominator sign):')
    for row in pat:
        print('   ', ' '.join(row))
