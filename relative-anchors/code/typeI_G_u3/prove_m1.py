"""Proof of G > 0 for u = |U| = 3, m = |T| = 1 (k = 4), all positive weights.

Input: build_m1.pkl (Bareiss leading minors of Nt = den G_z, each written as den^{e_j} c_j with c_j homogeneous in
(s0, s1, s2, s3)).  Since minor_j(G_z) = m_j / den^j = den^{e_j - j} c_j, G_z > 0 on the open orthant iff, for every
j, sign(den)^{e_j - j} sign(c_j) = +1 there (Sylvester's criterion).  Steps:
  (1) den: every irreducible factor has positive sign on the open orthant (all coefficients > 0, or a certificate);
  (2) each core c_j > 0 on the open orthant: either all its coefficients are positive, or, in the chart
      s0 = 1 - mu, s1 = mu a, s2 = mu (1 - a) b, s3 = mu (1 - a)(1 - b) ((mu, a, b) in (0,1)^3 covers every ray of the
      open orthant), after removing the largest powers of mu, a, b, 1 - a, 1 - b, 1 - mu, all tensor Bernstein
      coefficients on [0,1]^3 are positive (exact integers).
Every step is exact."""
import pickle, sys, time, flint
from bern import from_dict, certify

t00 = time.time()
R = pickle.load(open('build_m1.pkl', 'rb'))
names = R['names']
ctxZ = flint.fmpz_mpoly_ctx.get(tuple(names), 'lex')
ctx3 = flint.fmpz_mpoly_ctx.get(('mu', 'a', 'b'), 'lex')
mu, a, b = ctx3.gens()
SUB = (1 - mu, mu * a, mu * (1 - a) * b, mu * (1 - a) * (1 - b))
ok = True


def positive_orthant(p, name, max_boxes=4000):
    """True if p > 0 on the open orthant: all coefficients positive, or chart + strip + Bernstein."""
    coeffs = list(p.to_dict().values())
    if all(c > 0 for c in coeffs):
        print(f'  {name}: all {len(coeffs)} coefficients positive', flush=True)
        return True
    if all(c < 0 for c in coeffs):
        print(f'  {name}: all coefficients NEGATIVE', flush=True)
        return False
    q = p.compose(*SUB, ctx=ctx3)
    removed = []
    for nm, f in (('mu', mu), ('a', a), ('b', b), ('1-a', 1 - a), ('1-b', 1 - b), ('1-mu', 1 - mu)):
        e = 0
        while True:
            qq, r = divmod(q, f)
            if r != 0:
                break
            q = qq
            e += 1
        if e:
            removed.append((nm, e))
    dct = {tuple(int(x) for x in m): int(c) for m, c in q.to_dict().items()}
    degs = tuple(int(x) for x in q.degrees())
    t0 = time.time()
    good, st = certify(from_dict(dct, degs), max_boxes=max_boxes)
    print(f'  {name}: mixed signs; chart, removed {removed}, degrees (mu, a, b) = {degs}, {len(dct)} terms: '
          f'{"Bernstein > 0 on [0,1]^3" if good else "FAILED"} {st} ({time.time() - t0:.1f} s)', flush=True)
    return good


print('(1) the denominator den = cst * prod f_i^E_i', flush=True)
sgn_f = []
for k, (f, e) in enumerate(R['den_factors']):
    p = ctxZ.from_dict(f)
    if positive_orthant(p, f'den factor f_{k} (degree {p.total_degree()}) ^ {e}'):
        sgn_f.append(1)
    elif positive_orthant(-p, f'minus den factor f_{k}'):
        sgn_f.append(-1)
    else:
        sgn_f.append(0)
        ok = False
sgn_cst = 1 if R['den_const'] > 0 else -1
E = [e for f, e in R['den_factors']]
print(f'  signs of the factors on the open orthant: {sgn_f}; constant sign {sgn_cst}', flush=True)
print('(2) the cores of the leading minors: minor_j(G_z) = c_j prod f_i^(x_ij - j E_i) / cst^j', flush=True)
for j, (c, ex) in enumerate(R['cores']):
    jj = j + 1
    need = sgn_cst ** jj
    for s_, x, Ei in zip(sgn_f, ex, E):
        need *= s_ ** ((x - jj * Ei) % 2)
    p = ctxZ.from_dict(c)
    good = positive_orthant(need * p, f'core c_{jj} (exponents {ex}, degree {p.total_degree()}), required sign {need}')
    ok &= good
print(f'total {time.time() - t00:.1f} s')
print('PASS' if ok else 'FAIL')
