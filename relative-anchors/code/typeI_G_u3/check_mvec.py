"""Exact check (m = 1, u = 1, 2, 3): with v = (4 s_0 (s_0 + S), 3 (2 s_0 + S), ..., 3 (2 s_0 + S)) (s_0 = s_t; a positive
multiple of (1, beta, ..., beta), beta = 3 (2 s_t + S) / (4 s_t (s_t + S))), every entry of den * G_z v is a polynomial
with only positive coefficients.  Together with the Z-pattern (check_zpattern.py) this writes D_v G_z D_v as a weighted
graph Laplacian with subtraction-free edge weights plus the subtraction-free diagonal diag(v_i (G_z v)_i), so by the
all-minors matrix-forest theorem every principal minor of G_z is subtraction-free (Conjecture P for these cases, by a
second route), and G_z is a nonsingular M-matrix.  Floating-point evidence that the same v works for u <= 7: ex_beta.py,
ex_beta2.py (window for beta)."""
import pickle, flint

ok = True
for u in (1, 2, 3):
    D = pickle.load(open(f'closed_u{u}_m1.pkl', 'rb'))
    ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
    G = D['G']
    n = len(G)
    nums = [[ctx.from_dict(G[i][j][0]) for j in range(n)] for i in range(n)]
    den = ctx.from_dict(G[0][0][1])
    assert all(ctx.from_dict(G[i][j][1]) == den for i in range(n) for j in range(n))
    s = ctx.gens()
    s0, S = s[0], sum(s[1:])
    v = [4 * s0 * (s0 + S)] + [3 * (2 * s0 + S)] * u
    res = []
    for i in range(n):
        p = sum((nums[i][j] * v[j] for j in range(n)), 0 * den)
        c = list(p.to_dict().values())
        pos = all(x > 0 for x in c)
        res.append((pos, len(c), sum(1 for x in c if x < 0)))
        ok &= pos
    print(f'u = {u}: den * (G_z v)_i has only positive coefficients for every i: {all(r[0] for r in res)} '
          f'(terms, negative coefficients per row: {[(r[1], r[2]) for r in res]})', flush=True)
print('PASS' if ok else 'FAIL')
