"""Exact checks for Propositions 7 and 8 of NOTES.md, Section 5 (no anchor near the W-plane; anchors among single
excitations), n = 3..7.

Rows r_S (letter 1 at the two copies of S, 0 elsewhere), columns e_t (letter 1 at t).  p = sum_t w_t |e_t>.
G = M V_p, G[(a, y), x] = sum_t w_t M[(a, y), (e_t, x)], with M[(a,y),(b,x)] = prod_s (d(a_s,b_s) d(y_s,x_s) - 1/2 d(a_s,y_s) d(b_s,x_s)).
(i)   every row (r_S, e_t) of G is supported on the vacuum column x = 0...0, so on these states
      R_p^Gamma[(a,y'),(b,y)] = R_p[(b,y'),(a,y)] = G[(b,y'),0] G[(a,y),0] / D[0,0]   (D_p is block diagonal, Lemma 1);
(ii)  W_n (all w_t = 1): the Gram matrix of Y on (a_n, b_n) is 2^{1-n}(n-1)(n-2) [[0,-1],[-1,2(n+1)]];
(iii) rational weights: <a(w), Y a(w)> = 0, |a(w)|^2 = (n-1)|w|^2, <a(w), Y (r_ij (x) e_k)> = -2 c_n w_i w_j w_k / |w|^2,
      <b, Y b> = 4 c_n for b = r_ij (x) e_k, and b is orthogonal to a(w)  (c_n = 2^{1-n}).
All arithmetic in Fractions."""
import itertools
from fractions import Fraction as Fr

H = Fr(1, 2)

def m1(a, y, b, x): return (1 if (a == b and y == x) else 0) - (H if (a == y and b == x) else 0)

def Mel(a, y, b, x):
    v = Fr(1)
    for s in range(len(a)):
        v *= m1(a[s], y[s], b[s], x[s])
        if v == 0: return v
    return v

def Grow(a, y, w):
    """Nonzero entries {x: G[(a,y),x]} by per-copy enumeration (x_s in {y_s, b_s})."""
    n = len(a); out = {}
    for t in range(n):
        b = tuple(1 if s == t else 0 for s in range(n))
        for x in itertools.product(*[sorted({y[s], b[s]}) for s in range(n)]):
            v = Mel(a, y, b, x)
            if v: out[x] = out.get(x, 0) + w[t] * v
    return {x: v for x, v in out.items() if v != 0}

def run(n, w):
    vac = (0,) * n; cn = Fr(2) ** (1 - n)
    rS = lambda S: tuple(1 if s in S else 0 for s in range(n)); e = lambda t: tuple(1 if s == t else 0 for s in range(n))
    states = [(rS(S), e(t)) for S in itertools.combinations(range(n), 2) for t in range(n)]
    G = {st: Grow(st[0], st[1], w) for st in states}
    supp_ok = all(set(g) <= {vac} for g in G.values())
    D00 = sum(w[t] * w[u] * Mel(e(t), vac, e(u), vac) for t in range(n) for u in range(n))
    g0 = {st: G[st].get(vac, Fr(0)) for st in states}
    def Y(st1, st2):
        (a, yp), (b, y) = st1, st2
        return Mel(b, yp, a, y) - g0[(b, yp)] * g0[(a, y)] / D00
    def form(u, v): return sum(cu * cv * Y(s1, s2) for s1, cu in u.items() for s2, cv in v.items())
    return supp_ok, D00, states, Y, form, cn, rS, e

ok = True
for n in range(3, 8):
    w = [Fr(1)] * n
    supp_ok, D00, states, Y, form, cn, rS, e = run(n, w)
    a = {}; b = {}
    for S in itertools.combinations(range(n), 2):
        for t in range(n):
            (a if t in S else b)[(rS(S), e(t))] = Fr(1)
    Gm = [[form(a, a), form(a, b)], [form(b, a), form(b, b)]]
    pred = [[0, -cn * (n - 1) * (n - 2)], [-cn * (n - 1) * (n - 2), cn * (n - 1) * (n - 2) * 2 * (n + 1)]]
    good = supp_ok and D00 == n * cn and Gm == pred
    print(f"n={n} W: G supported on the vacuum column: {supp_ok}; D[0,0] = {D00} (= n c_n: {D00 == n * cn}); "
          f"Gram on (a_n, b_n) = {[[str(x) for x in r] for r in Gm]}: {Gm == pred}")
    ok &= good
    # weighted witness at two rational weight vectors
    for wv in ([Fr(i + 1, i + 2) for i in range(n)], [Fr(3), Fr(1, 2)] + [Fr(2 * i + 1, 7) for i in range(n - 2)]):
        supp_ok, D00, states, Y, form, cn, rS, e = run(n, wv)
        W2 = sum(x * x for x in wv)
        aw = {}
        for S in itertools.combinations(range(n), 2):
            for t in S:
                other = [x for x in S if x != t][0]; aw[(rS(S), e(t))] = wv[other]
        good = supp_ok and form(aw, aw) == 0 and sum(c * c for c in aw.values()) == (n - 1) * W2
        for i, j, k in itertools.permutations(range(n), 3):
            if i < j:
                bb = {(rS((i, j)), e(k)): Fr(1)}
                good &= form(aw, bb) == -2 * cn * wv[i] * wv[j] * wv[k] / W2 and form(bb, bb) == 4 * cn and (rS((i, j)), e(k)) not in aw
        print(f"n={n} weights {[str(x) for x in wv[:4]]}...: support on vacuum {supp_ok}; witness identities exact: {good}")
        ok &= good
print("PASS" if ok else "FAIL")
