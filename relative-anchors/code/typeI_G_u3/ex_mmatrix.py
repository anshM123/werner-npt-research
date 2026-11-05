"""Exploration (exact): for m = 1, is G_z a symmetric Z-matrix (off-diagonal entries <= 0) that is diagonally dominant,
i.e. a nonsingular M-matrix?  Prints the coefficient sign pattern of every entry and of the row sums G_z 1 (each row sum
written over the common denominator)."""
import pickle, sys, flint
for (u, m) in [(1, 1), (2, 1), (3, 1)]:
    D = pickle.load(open(f'closed_u{u}_m{m}.pkl', 'rb'))
    ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
    G = D['G']
    n = len(G)
    nums = [[ctx.from_dict(G[i][j][0]) for j in range(n)] for i in range(n)]
    dens = [[ctx.from_dict(G[i][j][1]) for j in range(n)] for i in range(n)]
    den = dens[0][0]
    same = all(d == den for row in dens for d in row)
    def sg(p):
        c = list(p.to_dict().values())
        return '+' if all(x > 0 for x in c) else ('-' if all(x < 0 for x in c) else ('0' if not c else '?'))
    print(f'(u, m) = ({u}, {m}): common denominator: {same}, den sign {sg(den)}')
    for i in range(n):
        print('   entries:', ' '.join(sg(nums[i][j]) for j in range(n)), '  row sum numerator sign:', sg(sum((nums[i][j] for j in range(n)), 0 * den)) if same else 'n/a')
