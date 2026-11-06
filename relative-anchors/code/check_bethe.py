"""Exact check of the Bethe ansatz for E_c on strings with the letters 0 and 1 (Lemma BA of the README, Section 2).
With B(x) = sum_t w_t sigma_t^+ / (s_t - x) and |x> = B(x_1) ... B(x_c)|0>, the off-shell identity

    E |x> = (2|s| - 2 sum_a x_a) |x> + sum_a (Lambda_a / 2 - 2) J^+ |x^(a)>,
    Lambda_a = sum_t s_t/(s_t - x_a) + 2 sum_{b != a} x_b/(x_a - x_b),   |x^(a)> = prod_{b != a} B(x_b)|0>,

is checked as an identity of rational functions in the symbols x_1, ..., x_c (sympy), with the weights w_t symbolic as
well (s_t = w_t^2), for every 1 <= c <= k <= 4, and with rational weights for k = 5, c <= 3.  E = C + J^+J^-/2 with
C = 2 sum_t s_t |0><0|_t is built on the strings directly.  python check_bethe.py"""
import itertools
import sympy as sp


def run(k, c, symbolic):
    if symbolic:
        w = sp.symbols(f'w0:{k}', positive=True)
    else:
        w = [sp.Rational(n, d) for n, d in [(1, 1), (2, 3), (3, 2), (1, 2), (5, 4)][:k]]
    s = [x ** 2 for x in w]
    xs = sp.symbols(f'x0:{c}')
    strings = {cc: [frozenset(P) for P in itertools.combinations(range(k), cc)] for cc in range(c + 1)}

    def Bvec(roots):
        """prod_a B(x_a)|0> as a dict over c-subsets"""
        v = {frozenset(): sp.Integer(1)}
        for x in roots:
            nv = {}
            for y, co in v.items():
                for t in range(k):
                    if t not in y:
                        z = y | {t}
                        nv[z] = nv.get(z, 0) + co * w[t] / (s[t] - x)
            v = nv
        return v

    def Jp(v):
        out = {}
        for y, co in v.items():
            for t in range(k):
                if t not in y:
                    z = y | {t}
                    out[z] = out.get(z, 0) + co * w[t]
        return out

    def Jm(v):
        out = {}
        for y, co in v.items():
            for t in y:
                z = y - {t}
                out[z] = out.get(z, 0) + co * w[t]
        return out

    def Eop(v):
        out = {}
        for y, co in v.items():
            out[y] = out.get(y, 0) + co * 2 * sum(s[t] for t in range(k) if t not in y)
        for y, co in Jp(Jm(v)).items():
            out[y] = out.get(y, 0) + co / 2
        return out

    X = Bvec(xs)
    lhs = Eop(X)
    S = sum(s)
    rhs = {y: (2 * S - 2 * sum(xs)) * co for y, co in X.items()}
    for a in range(c):
        Lam = sum(s[t] / (s[t] - xs[a]) for t in range(k)) + 2 * sum(xs[b] / (xs[a] - xs[b]) for b in range(c) if b != a)
        Xa = Bvec([xs[b] for b in range(c) if b != a])
        for y, co in Jp(Xa).items():
            rhs[y] = rhs.get(y, 0) + (Lam / 2 - 2) * co
    bad = 0
    for y in strings[c]:
        if sp.simplify(sp.together(lhs.get(y, 0) - rhs.get(y, 0))) != 0:
            bad += 1
    return bad


ok = True
for k in range(1, 5):
    for c in range(1, k + 1):
        b = run(k, c, True)
        print(f'k={k} c={c} (symbolic weights): {b} strings where the identity fails', flush=True)
        ok &= b == 0
for c in range(1, 4):
    b = run(5, c, False)
    print(f'k=5 c={c} (rational weights): {b} strings where the identity fails', flush=True)
    ok &= b == 0
print('PASS' if ok else 'FAIL')
