"""Proposition E2 (case II letter-2 blocks, every k, j at equal weights): closed form
   Q / |psi|^2 = 1/2 - (1/4) [ a (b-1)^2 / ((a+b) e_plus) + b (a+1)^2 / ((a+b) e_minus) ],
   a = q + 1, b = m - q, e_plus = (4b + a + ab)/2 + kap, e_minus = b(a + 3)/2 + kap, kap = 2i + tau1,
   q = j - 1 - tau1 - i (Dicke count), m = k - 1 - tau1 - 2i.
Checks the closed form against eqlib2 (explicit vectors, Lemma B) exactly for k = 2..9, every j, tau1, i."""
from fractions import Fraction
from eqlib2 import caseII_vectors, form

ok = True
n = 0
for k in range(2, 10):
    for j in range(0, k + 1):
        for tau1 in range(0, k):
            for i in range(0, (k - 1 - tau1) // 2 + 1):
                vecs = caseII_vectors(k, j, tau1, i)
                if not vecs:
                    continue
                G, Q = form(k, j, vecs)
                m = k - 1 - tau1 - 2 * i
                q = j - 1 - tau1 - i
                a, b, kap = q + 1, m - q, 2 * i + tau1
                ep = Fraction(4 * b + a + a * b, 2) + kap
                em = Fraction(b * (a + 3), 2) + kap
                val = Fraction(1, 2) - Fraction(1, 4) * (Fraction(a * (b - 1) ** 2, (a + b)) / ep
                                                       + (Fraction(b * (a + 1) ** 2, (a + b)) / em if b > 0 else 0))
                got = Fraction(int(Q[0][0].p), int(Q[0][0].q)) / Fraction(int(G[0][0].p), int(G[0][0].q))
                n += 1
                if got != val:
                    ok = False
                    print('MISMATCH', k, j, tau1, i, got, val)
print(n, 'case II forms compared;', 'PASS' if ok else 'FAIL')
