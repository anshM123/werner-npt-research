"""EXACT check (rational point, rational w) of Lemma G (relative-anchors/README.md): on G0_{j+1} (row = vacuum, column in H_{j+1}, all
letters 0,1,2), UNRESCALED,
    Y00 = 2^{-(k-j-1)} (I - 1/2 J^- E^{-1} J^+),   E = C + 1/2 J^+ J^-  on H_{j+2},
    C = diag( sum_t s_t (2[y_t=0] + [y_t=2]) ),
and, when C is invertible on H_{j+2} (j + 2 <= k - 1),  Y00^{-1} = 2^{k-j-1} (I + 1/2 J^- C^{-1} J^+)  (Woodbury).
Also the letter-1 blocks of (R_inf)_0 restricted to G1_0 and of (R_inf)_k (RESCALED):
    G1_0 letter-1 block = 2^{-(k-1)} (diag(s) - s s^T/|s|),   (R_inf)_k letter-1 block = 1/2 prod(s) (diag(s) - s s^T/|s|).
k = 2..4.  Prints PASS/FAIL."""
import itertools, random
import flint
from ylib import Builder, point_ring, Nv, inv
Q_ = flint.fmpq
random.seed(3)
ok = True
for k in (2, 3, 4):
    w = [Q_(random.randint(1, 9), random.randint(1, 9)) for _ in range(k)]
    s = [x * x for x in w]; S = sum(s, Q_(0))
    Bd = Builder(point_ring(s))
    def mu(X):
        v = Q_(1)
        for t, n in enumerate(Nv(*X)):
            v *= w[t] ** n
        return v
    q = lambda y: sum(x != 0 for x in y)
    allc = list(itertools.product(range(3), repeat=k))
    vac = (0,) * k
    for j in range(-1, k):
        H1 = [y for y in allc if q(y) == j + 1]; H2 = [y for y in allc if q(y) == j + 2]
        st = [(vac, y) for y in H1]
        Y00 = [[Bd.Y(X1, X2) / (mu(X1) * mu(X2)) for X2 in st] for X1 in st]
        # J^+ : H1 -> H2, letter 1 at an empty copy
        Jp = [[sum((w[t] for t in range(k) if y[t] == 0 and y2 == tuple(1 if u == t else y[u] for u in range(k))), Q_(0))
               for y in H1] for y2 in H2]
        Cd = [sum((s[t] * (2 * (y[t] == 0) + (y[t] == 2)) for t in range(k)), Q_(0)) for y in H2]
        n1, n2 = len(H1), len(H2)
        JpJm = [[sum((Jp[a][c] * Jp[b][c] for c in range(n1)), Q_(0)) for b in range(n2)] for a in range(n2)]
        E = [[(Cd[a] if a == b else Q_(0)) + JpJm[a][b] / 2 for b in range(n2)] for a in range(n2)]
        Ei = inv(E, Q_(0), Q_(1)) if n2 else []
        pref = Q_(1, 2 ** (k - j - 1))
        pred = [[pref * ((Q_(1) if a == b else Q_(0)) - sum((Jp[c][a] * Ei[c][d] * Jp[d][b] for c in range(n2) for d in range(n2)), Q_(0)) / 2)
                 for b in range(n1)] for a in range(n1)]
        g1 = all(Y00[a][b] == pred[a][b] for a in range(n1) for b in range(n1))
        g2 = True
        if n2 and all(c != 0 for c in Cd):
            Yi = inv(Y00, Q_(0), Q_(1))
            pi = [[(2 ** (k - j - 1)) * ((Q_(1) if a == b else Q_(0)) + sum((Jp[c][a] * Jp[c][b] / Cd[c] for c in range(n2)), Q_(0)) / 2)
                   for b in range(n1)] for a in range(n1)]
            g2 = all(Yi[a][b] == pi[a][b] for a in range(n1) for b in range(n1))
            tag = 'Woodbury inverse checked'
        else:
            tag = 'C singular or H_{j+2} empty: inverse formula not applicable'
        ok &= g1 and g2
        print(f'k={k} j={j}: Y00 = 2^-(k-j-1)(I - J^-E^-1J^+/2): {g1}; {tag}: {g2}')
    # letter-1 blocks
    B1 = Builder(point_ring(s), letters=(0, 1))
    e = [tuple(1 if r == t else 0 for r in range(k)) for t in range(k)]
    st = [(e[r], vac) for r in range(k)]
    M0 = [[B1.Y(X1, X2) for X2 in st] for X1 in st]
    g3 = all(M0[a][b] == Q_(1, 2 ** (k - 1)) * ((s[a] if a == b else Q_(0)) - s[a] * s[b] / S) for a in range(k) for b in range(k))
    one = (1,) * k
    st = [(e[r], one) for r in range(k)]
    Mk = [[B1.Y(X1, X2) for X2 in st] for X1 in st]
    P = Q_(1)
    for x in s: P *= x
    g4 = all(Mk[a][b] == P / 2 * ((s[a] if a == b else Q_(0)) - s[a] * s[b] / S) for a in range(k) for b in range(k))
    ok &= g3 and g4
    print(f'k={k}: G1_0 letter-1 Laplacian form {g3}; (R_inf)_k letter-1 Laplacian form {g4}')
print('PASS' if ok else 'FAIL')
