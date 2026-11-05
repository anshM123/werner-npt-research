"""Test 2 (exact negative statement): no representation with fixed vectors.

If M(s) = sum_e c_e(s) x_e x_e^T with c_e >= 0 on the open orthant and fixed real vectors x_e (any finite family),
then |M_ab| <= K_ab M_aa with K_ab = max_e |x_e(b) / x_e(a)| (over e with x_e(a) != 0).  So M_ab / M_aa must be bounded
on the open orthant, for every ordered pair (a, b).  With M_ab = n_ab / d_ab and M_aa = n_aa / d_aa, where d_ab and n_aa
have positive coefficients, the ratio P/Q (P = n_ab d_aa, Q = d_ab n_aa) is bounded iff conv(supp P) lies in
conv(supp Q).  An exact certificate of unboundedness is an integer direction l with
max_{m in supp P} l.m > max_{m in supp Q} l.m: along s = (t^{l_0}, ..., t^{l_r}), t -> infinity, the ratio grows like
t^{gap} (the maximum over supp P is attained at a vertex, whose coefficient cannot cancel, after a small perturbation of l
that we also check: we require the maximum over supp P to be attained at a single monomial).
For each block the script reports the pairs (a, b) with an exact certificate.  Pairs without a certificate among the
directions in {-3..3}^r are tested for boundedness by linear programming (floating point, reported separately; skipped
with the argument nolp).
Width certificate: a rescaling of the coordinates by monomials, x_a -> s^{delta_a} x_a, multiplies M_ab / M_aa by the
monomial s^{delta_b - delta_a}, i.e. translates conv(supp P).  If conv(supp P) is wider than conv(supp Q) in some
integer direction l (max - min of l.m over supp P larger than over supp Q), no translate fits, so the obstruction
survives every monomial rescaling (exact).
Usage: python newton_obstr.py k [nolp]"""
import sys, itertools, time
import numpy as np
from scipy.optimize import linprog
from psdtools import load_blocks

k = int(sys.argv[1])
SKIPLP = len(sys.argv) > 2 and sys.argv[2] == 'nolp'
t0 = time.time()
tot_blocks = 0
obstructed = 0
wobstructed = 0
for blk in load_blocks(k):
    n = len(blk['Mg'])
    Mg = blk['Mg']
    nv = len(blk['names'])
    dirs = [l for l in itertools.product(range(-3, 4), repeat=nv) if any(l)]
    dirs.sort(key=lambda l: sum(abs(x) for x in l))
    certs = []
    lp_unb = []
    for a in range(n):
        na, da = Mg[a][a]
        assert all(c > 0 for c in na.coeffs()), 'diagonal numerator not positive'
        for b in range(n):
            if b == a or Mg[a][b][0].is_zero():
                continue
            nab, dab = Mg[a][b]
            P = nab * da
            Qp = dab * na
            if not all(c > 0 for c in Qp.coeffs()):
                Qp = -Qp
            assert all(c > 0 for c in Qp.coeffs())
            sP = list(P.to_dict().keys())
            sQ = list(Qp.to_dict().keys())
            found = None
            for l in dirs:
                vP = [sum(x * y for x, y in zip(l, m)) for m in sP]
                vQ = max(sum(x * y for x, y in zip(l, m)) for m in sQ)
                mp = max(vP)
                if mp > vQ and vP.count(mp) == 1:
                    found = (l, mp - vQ)
                    break
            if found:
                # width certificate: Newt(P) wider than Newt(Q) in some direction l, so no translation of Newt(P)
                # (i.e. no rescaling of the coordinates by monomials) fits it into Newt(Q)
                wfound = None
                for l in dirs:
                    vP = [sum(x * y for x, y in zip(l, m)) for m in sP]
                    vQ = [sum(x * y for x, y in zip(l, m)) for m in sQ]
                    if max(vP) - min(vP) > max(vQ) - min(vQ):
                        wfound = l
                        break
                certs.append((a, b, found, wfound))
            elif not SKIPLP:
                # LP: is every monomial of P in conv(supp Q)?
                A_eq = np.array([[float(m[i]) for m in sQ] for i in range(nv)] + [[1.0] * len(sQ)])
                unb = False
                for m in sP:
                    res = linprog(np.zeros(len(sQ)), A_eq=A_eq, b_eq=np.array([float(x) for x in m] + [1.0]),
                                  bounds=[(0, None)] * len(sQ), method='highs')
                    if res.status != 0:
                        unb = True
                        break
                if unb:
                    lp_unb.append((a, b))
    tot_blocks += 1
    if certs:
        obstructed += 1
    labs = blk['labels']
    wcerts = [c for c in certs if c[3] is not None]
    if wcerts:
        wobstructed += 1
    msg = (f'{blk["name"]} ({blk["kind"]}), size {n}: ordered pairs with M_ab/M_aa unbounded, exact certificate: '
           f'{len(certs)} (of which {len(wcerts)} also with a width certificate); LP-only: '
           f'{"skipped" if SKIPLP else len(lp_unb)}')
    if certs:
        a, b, (l, gap), wl = certs[0]
        msg += f'; e.g. M[{labs[a]},{labs[b]}]/M[{labs[a]},{labs[a]}] ~ t^{gap} along s = t^{l}'
    if wcerts:
        a, b, _, wl = wcerts[0]
        msg += f'; width certificate for M[{labs[a]},{labs[b]}]/M[{labs[a]},{labs[a]}] in direction {wl}'
    print(msg, flush=True)
    print(f'   => fixed-vector representation {"IMPOSSIBLE (exact)" if certs else "not excluded by this test"}; '
          f'with any monomial rescaling of the coordinates: {"IMPOSSIBLE (exact)" if wcerts else "not excluded"}')
print(f'RESULT k = {k}: {obstructed} of {tot_blocks} blocks have an exact obstruction, {wobstructed} of them also '
      f'against every monomial rescaling  [{time.time() - t0:.1f}s]')
print('PASS')
