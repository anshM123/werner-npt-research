"""Negative controls for the determinant routes (NOTES.md, Section 9): B -> (1 + 1/1000) B in the calculus.

For each family (sym mode): the form with B scaled (tfamlib bscale), the split (block-diagonality still checked, it does
not depend on the scaling), and, with the same gauge as the true form:
  (a) <kappa, Q kappa> for each known kernel vector at the minimal member and sigma = 1 (exact; negative means the
      scaled form is not >= 0 there, so no certificate can succeed);
  (b) the determinant of the gauge-fixed first part at the minimal member as a polynomial in sigma (exact): its
      numerator after clearing the (one-signed) denominator, the number of negative coefficients, and its sign at
      sigma = 10^-6, 1, 10^6.
The control passes (the routes would fail, as they must) when some family has (a) negative or (b) a negative value.
python negdet.py FAM [FAM ...]     prints FAIL (expected) when the scaled forms are caught, PASS otherwise"""
import sys, time
from fractions import Fraction as Fr
import flint
import tfamlib as T
import tprove as P
import tsplit
import tsbuild


def control(case, fam, b=Fr(1001, 1000), log=print):
    C = T.FamCalc(case, 'sym', *fam, bscale=b)
    names, Q = T.reduced(C)
    vecs, two, three, split = tsplit.basis(case, names, C.q0, C.m, C.i, C.one)
    M = tsplit.congruence(names, Q, vecs, two + three, C.zero)
    blockdiag = all(M[p, q] == 0 for p in two for q in three)
    kv = [(lab, tsplit.to_new(v, split)) for lab, v in T.kernel_vectors(C, names)]
    mp = C.minpoint()
    pt = [mp[nm] for nm in C.syms if nm != 's']
    kqk = []
    for lab, w in kv:
        s = C.zero
        for a, ca in w.items():
            for c, cc in w.items():
                s = s + ca * cc * M[a, c]
        kqk.append(P.evalF(s, pt + [1]))
    # the same gauge as the true form
    Ct = T.FamCalc(case, 'sym', *fam)
    namest, Qt = T.reduced(Ct)
    vt, twot, threet, splitt = tsplit.basis(case, namest, Ct.q0, Ct.m, Ct.i, Ct.one)
    gauge = []
    if kv:
        import itertools
        pref = [nm for nm in tsbuild.GPREF[case] if nm in two] + [nm for nm in two if nm not in tsbuild.GPREF[case]]
        for S in itertools.combinations(pref, len(kv)):
            sub = [[C.K(w.get(nm, 0)) for nm in S] for _, w in kv]
            d = P.det_small(sub, C)
            if d != 0 and d.n.is_constant() and d.d.is_constant():
                gauge = list(S)
                break
    keep = [nm for nm in two if nm not in gauge]
    # determinant at the minimal member, sigma free
    K1 = flint.fmpq_poly
    vals = {nm: flint.fmpq(mp[nm]) for nm in C.syms if nm != 's'}

    def at(x):
        if x == 0:
            return (K1([0]), K1([1]))
        args = []
        sp = flint.fmpq_poly([0, 1])
        nn = x.n.to_dict()
        dd = x.d.to_dict()
        idx = [C.syms.index(nm) for nm in C.syms if nm != 's']
        isg = C.syms.index('s')

        def ev(dct):
            out = K1([0])
            for e, c in dct.items():
                term = flint.fmpq(c)
                for k_ in idx:
                    term *= vals[C.syms[k_]] ** e[k_]
                out += K1([0] * e[isg] + [term])
            return out
        return (ev(nn), ev(dd))
    n = len(keep)
    A = [[at(M[a, c]) for c in keep] for a in keep]
    # fraction-free determinant over Q[sigma]: Bareiss on the cleared matrix
    den = K1([1])
    for row in A:
        for (p, q) in row:
            den = den * q // den.gcd(q) if q != 0 else den
    Mx = [[p * (den // q) for (p, q) in row] for row in A]
    prev = K1([1])
    sign = 1
    for k_ in range(n - 1):
        if Mx[k_][k_] == 0:
            r = next((r for r in range(k_ + 1, n) if Mx[r][k_] != 0), None)
            if r is None:
                return {'det': 0}
            Mx[k_], Mx[r] = Mx[r], Mx[k_]
            sign = -sign
        for i_ in range(k_ + 1, n):
            for j_ in range(k_ + 1, n):
                Mx[i_][j_] = (Mx[i_][j_] * Mx[k_][k_] - Mx[i_][k_] * Mx[k_][j_]) // prev
        prev = Mx[k_][k_]
    detp = Mx[n - 1][n - 1] * sign
    co = [detp[i] for i in range(detp.degree() + 1)]
    nneg = sum(1 for c in co if c < 0)
    npos = sum(1 for c in co if c > 0)
    dsgn = 1 if den(flint.fmpq(1)) > 0 else -1
    vals_s = {}
    for sv in (flint.fmpq(1, 10 ** 6), flint.fmpq(1), flint.fmpq(10 ** 6)):
        v = detp(sv) * dsgn ** n
        vals_s[str(sv)] = 1 if v > 0 else (-1 if v < 0 else 0)
    return {'blockdiag': blockdiag, 'kQk': kqk, 'gauge': gauge, 'det_coeffs': (npos, nneg), 'det_signs': vals_s,
            'n': n}


if __name__ == '__main__':
    T0 = time.time()
    caught = 0
    for spec in sys.argv[1:]:
        case, mode, fam = P.parse_fam(spec)
        r = control(case, fam)
        neg = any(x < 0 for x in r.get('kQk', [])) or any(v < 0 for v in r.get('det_signs', {}).values())
        caught += neg
        print(f'{P.fam_label(case, mode, fam)} B scaled by 1001/1000: split block-diagonal {r.get("blockdiag")}; '
              f'<kappa, Q kappa> at the minimal member, sigma = 1: {[float(x) for x in r.get("kQk", [])]}; '
              f'det of the gauge-fixed first part ({r.get("n")} x {r.get("n")}, gauge {r.get("gauge")}) at the '
              f'minimal member as a polynomial in sigma: {r.get("det_coeffs")} (positive, negative) coefficients, '
              f'signs at sigma = 1e-6, 1, 1e6: {r.get("det_signs")}; caught: {neg} [{time.time()-T0:.0f}s]',
              flush=True)
    print(f'negative control: {caught} of {len(sys.argv) - 1} scaled families caught (not >= 0, or a negative '
          f'determinant of the gauge-fixed part)')
    print('FAIL (expected: the scaled forms are caught)' if caught > 0 else 'PASS (nothing caught)')
