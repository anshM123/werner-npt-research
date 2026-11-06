"""Conjecture R at equal U-weights, symbolically in (sigma, s, u) (python-flint rational functions).

With s_h = s for all h, Psi commutes with the permutations of U, so (NOTES.md, Section 9)
  Psi > 0 on {a orthogonal to w}  <=>  D(s) > 0 (standard part, u >= 2)  and  Psi_inv > 0,
where Psi_inv is Psi on the invariant vectors (ah_h = 0 by the constraint, ch_h = ch, zh_h = zh) and rho (a 3 x 3
matrix).  All entries are rational functions of (sigma, s, u); S = u s, M_k = u phi_k(s), IG = 1/(2 + u s/d),
IGT = 1/(2 + sigma/d0_t + u s/d0).  The script computes the leading principal minors of D(s) and Psi_inv and tests,
after the substitution u = u0 + x, whether every numerator and every irreducible denominator factor has coefficients
of one sign as a polynomial in (sigma, s, x).  Output: equal_R.pkl (the six minors)."""
import sys, time, pickle
from fractions import Fraction as Fr
import sympy as sp
import flint
from flintfield import Field
from reduced_c import ReducedC as Reduced

s_, sig_, S_ = sp.symbols('s sigma S', positive=True)


def main(u0):
    t0 = time.time()
    R = Reduced()
    Kf = Field(['sigma', 's', 'x'])
    sg, sv, x = Kf.gens
    u = x + u0
    S = u * sv
    d = sv + S / 2 + sg
    d0 = sv + (sg + S) / 2
    dt0 = (3 * sg + S) / 2
    IGv = Kf(1) / (2 + u * sv / d)
    IGTv = Kf(1) / (2 + sg / dt0 + u * sv / d0)

    def pw(v, e):
        r = Kf(1)
        for _ in range(e):
            r = r * v
        return r

    def ev_terms(terms, vals):
        tot = Kf(0)
        for m, c in terms:
            t = Kf(c if isinstance(c, flint.fmpq) else flint.fmpq(c.numerator, c.denominator))
            for v, e in zip(vals, m):
                if e:
                    t = t * pw(v, e)
            tot = tot + t
        return tot

    def basis_val(tn_td):
        return ev_terms(tn_td[0], (sv, sg, S)) / ev_terms(tn_td[1], (sv, sg, S))
    Mvals = [u * basis_val(t) for t in R.Mbasis]
    cv = [IGv, IGTv, u] + Mvals
    print(f'constants built ({time.time() - t0:.1f} s)', flush=True)

    cache = {}

    def cmon(ex):
        if ex not in cache:
            t = Kf(1)
            for v, e in zip(cv, ex):
                if e:
                    t = t * pw(v, e)
            cache[ex] = t
        return cache[ex]
    # D(s)
    names = ['ah', 'ch', 'zh']
    D3 = [[Kf(0)] * 3 for _ in range(3)]
    for a, b, ex, nt, dt in R.loc:
        val = cmon(ex) * ev_terms(nt, (sv, sg, S)) / ev_terms(dt, (sv, sg, S))
        ia, ib = names.index(a), names.index(b)
        if ia == ib:
            D3[ia][ia] = D3[ia][ia] + val
        else:
            D3[ia][ib] = D3[ia][ib] + val / 2
            D3[ib][ia] = D3[ib][ia] + val / 2
    print(f'D(s) built ({time.time() - t0:.1f} s)', flush=True)
    # Psi_inv on (ch, zh, rho): y_0 = rho; y_j = u phi_j(s) * (ch or zh) for F_j of class ch, zh; 0 for class ah
    L = [(2, Kf(1))]
    for xcls, t in R.Fbasis:
        if xcls == 'ah':
            L.append(None)
        else:
            L.append((names.index(xcls) - 1, u * basis_val(t)))
    P3 = [[Kf(0)] * 3 for _ in range(3)]
    for i in range(2):
        for j in range(2):
            P3[i][j] = u * D3[1 + i][1 + j]
    nterm = 0
    for i, j, ex, nt, dt in R.Cterms:
        if L[i] is None or L[j] is None:
            continue
        val = cmon(ex) * ev_terms(nt, (sg, S)) / ev_terms(dt, (sg, S))
        p, li = L[i]
        q, lj = L[j]
        contrib = val * li * lj
        if i == j:
            P3[p][q] = P3[p][q] + contrib
        else:
            P3[p][q] = P3[p][q] + contrib / 2
            P3[q][p] = P3[q][p] + contrib / 2
        nterm += 1
    print(f'Psi_inv built from {nterm} C-terms ({time.time() - t0:.1f} s)', flush=True)

    def minors(A):
        m1 = A[0][0]
        m2 = A[0][0] * A[1][1] - A[0][1] * A[1][0]
        m3 = (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
              + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
        return [m1, m2, m3]
    out = {}
    allok = True
    for name, A in (('D', D3), ('Psi_inv', P3)):
        for r, m in enumerate(minors(A), 1):
            num = m.n
            den = m.d
            cn = [c for c in num.coeffs()]
            sgn_n = 1 if all(c > 0 for c in cn) else (-1 if all(c < 0 for c in cn) else 0)
            fac = den.factor()
            sgn_d = []
            for f, e in fac[1]:
                cf = list(f.coeffs())
                sgn_d.append(1 if all(c > 0 for c in cf) else (-1 if all(c < 0 for c in cf) else 0))
            cst = fac[0]
            # overall sign of the denominator: constant times product of factor signs ** exponents
            dsign = (1 if cst > 0 else -1)
            for (f, e), sg_ in zip(fac[1], sgn_d):
                dsign *= sg_ ** e if sg_ != 0 else 0
            ok = sgn_n != 0 and dsign != 0 and sgn_n * dsign > 0
            allok &= ok
            print(f'{name}: leading minor {r}: numerator {len(cn)} terms, degree {num.total_degree()}, sign '
                  f'{sgn_n}; denominator factors {len(fac[1])} (signs {sgn_d}); positive certificate: {ok}',
                  flush=True)
            out[(name, r)] = (num.to_dict(), den.to_dict())
    pickle.dump(out, open(f'equal_R_u{u0}.pkl', 'wb'))
    print(f'total {time.time() - t0:.1f} s')
    print('PASS' if allok else 'FAIL')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
