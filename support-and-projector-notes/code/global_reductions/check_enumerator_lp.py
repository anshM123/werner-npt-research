"""Weight-enumerator identities for rank-two projectors and the linear-programming certificates (n = 4..9)."""
import itertools
from fractions import Fraction as Fr
from math import comb
import numpy as np
from qtools import q_dense, weyl, rand_c

ok = True
LAB = [(a, b) for a in range(3) for b in range(3)]


def enumerators(Pi, n):
    A = np.zeros(n + 1); B = np.zeros(n + 1)
    for labs in itertools.product(LAB, repeat=n):
        E = weyl(list(labs))
        w = sum(1 for l in labs if l != (0, 0))
        A[w] += abs(np.trace(Pi @ E)) ** 2 / 4
        B[w] += np.trace(Pi @ E @ Pi @ E.conj().T).real / 2
    return A, B


def poly_coeffs(A, n, f):
    """coefficients c_m of x^{n-m} y^m in sum_j A_j (f(x,y) substituted), f gives ((p,q),(r,s)) for
    (p x + q y)^{n-j} (r x + s y)^j."""
    (p, q), (rr, s) = f
    out = np.zeros(n + 1)
    for j in range(n + 1):
        for h in range(n - j + 1):          # choose y from the first factor h times
            for g in range(j + 1):          # choose y from the second factor g times
                out[h + g] += A[j] * comb(n - j, h) * p ** (n - j - h) * q ** h * comb(j, g) * rr ** (j - g) * s ** g
    return out


rng = np.random.default_rng(3)
for n in (2, 3):
    for trial in range(2):
        Qm, _ = np.linalg.qr(rand_c((3 ** n, 2), rng))
        Pi = Qm @ Qm.conj().T
        A, B = enumerators(Pi, n)
        qv = q_dense(Pi, n).real
        ok &= abs(qv - 4 / 3 ** n * sum(A[j] * (-0.5) ** (n - j) for j in range(n + 1))) < 1e-10
        ok &= abs(A.sum() - 3 ** n / 2) < 1e-9 and abs(A[0] - 1) < 1e-12
        Bm = 2 / 3 ** n * poly_coeffs(A, n, ((1, 8), (1, -1)))
        ok &= np.allclose(Bm, B)
        Sh = poly_coeffs(A, n, ((4, 2), (1, -1)))
        ok &= (Sh > -1e-9).all()
        # shadow polynomial = (6^n/4) Tr[(Pi (x) Pi) (x)_nu (x P+ + y P-)] at (x, y) = (0.3, 1.7)
        x, y = 0.3, 1.7
        lhs = 0.0
        for S in itertools.product((0, 1), repeat=n):
            # (x P+ + y P-) = (x+y)/2 I + (x-y)/2 F ; Tr[(Pi(x)Pi) F_S] = Tr[(Tr_{S^c} Pi)^2]
            T = Pi.reshape((3,) * (2 * n))
            for s in reversed([s for s in range(n) if not S[s]]):
                T = np.trace(T, axis1=s, axis2=s + T.ndim // 2)
            m = sum(S)
            Ms = T.reshape(3 ** m, 3 ** m) if m else T.reshape(1, 1)
            lhs += ((x + y) / 2) ** (n - m) * ((x - y) / 2) ** m * np.trace(Ms @ Ms).real
        rhs = 4 / 6 ** n * sum(Sh[mm] * x ** (n - mm) * y ** mm for mm in range(n + 1))
        ok &= abs(lhs - rhs) < 1e-9 * max(1, abs(lhs))
print('enumerator identities (n = 2, 3, random rank-two projectors):', ok)

# balanced three-copy code span{g, h}
g = np.zeros(27, complex); h = np.zeros(27, complex)
for a in range(3):
    g[9 * a + 3 * a + a] = 1 / np.sqrt(3)
for (a, b, c) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    h[9 * a + 3 * b + c] = 1 / np.sqrt(3)
Pi = np.outer(g, g.conj()) + np.outer(h, h.conj())
A, B = enumerators(Pi, 3)
ok &= abs(A[1]) < 1e-12 and abs(B[1]) < 1e-12 and abs(q_dense(Pi, 3).real - 1.5) < 1e-12
print('three-copy balanced code: A1 = B1 = 0, q_3 = 3/2:', abs(A[1]) < 1e-12, abs(B[1]) < 1e-12, q_dense(Pi, 3).real)


# --- exact linear-programming certificates ---
def enum_lp(n):
    d = n - 3
    empty = lambda: [Fr(0)] * (d + 1)
    A = [empty() for _ in range(n + 1)]
    M = Fr(3 ** n, 2); moment = Fr(8 * n, 9) * M
    A[0][0] = Fr(1)
    for j in range(2, n - 1):
        A[j][j - 1] = Fr(1)
    A[n - 1][0] = n * (M - 1) - moment
    A[n][0] = moment - (n - 1) * (M - 1)
    for j in range(2, n - 1):
        A[n - 1][j - 1] = Fr(j - n); A[n][j - 1] = Fr(n - 1 - j)
    comb_ = lambda rows: [sum(c * row[i] for c, row in rows) for i in range(d + 1)]
    cons = []
    for j in range(2, n + 1):
        cons.append(('A%d' % j, A[j]))
    for m in range(2, n + 1):
        rows = []
        for j in range(n + 1):
            z = sum(comb(n - j, hh) * 8 ** hh * comb(j, m - hh) * (-1) ** (m - hh) for hh in range(m + 1))
            rows.append((Fr(2 * z, 3 ** n), A[j]))
        Bm = comb_(rows)
        cons.append(('B%d-A%d' % (m, m), [Bm[i] - A[m][i] for i in range(d + 1)]))
    for m in range(n + 1):
        rows = []
        for j in range(n + 1):
            z = sum(comb(n - j, hh) * 4 ** (n - j - hh) * 2 ** hh * comb(j, m - hh) * (-1) ** (m - hh)
                    for hh in range(m + 1) if hh <= n - j and m - hh <= j)
            rows.append((Fr(z), A[j]))
        cons.append(('Shadow%d' % m, comb_(rows)))
    obj = comb_([(Fr((-1) ** (n - j), 2 ** (n - j)), A[j]) for j in range(n + 1)])
    return A, cons, obj


DATA = {
    4: ("297/16", ["1", "0", "0", "14", "51/2"], [0], ["9/4"], "11/12"),
    5: ("1215/32", ["1", "0", "0", "5", "105/2", "63"], [0, 4], ["0", "27/8"], "5/8"),
    6: ("6561/64", ["1", "0", "0", "7", "63/2", "153", "172"], [0, 5, 16], ["0", "81/16", "1/8"], "9/16"),
    7: ("2709693/21376", ["1", "0", "693/167", "3857/167", "7707/167", "0", "197701/334", "71340/167"],
        [3, 6, 8, 18], ["729/1336", "67959/5344", "891/2672", "51/2672"], "1239/5344"),
    8: ("14187069/60544", ["1", "0", "15233/946", "176621/3784", "331055/1892", "0", "0", "7106735/3784",
                           "2201615/1892"], [3, 4, 7, 9, 22],
        ["50301/30272", "19683/60544", "1400409/60544", "62937/60544", "6517/60544"], "6487/45408"),
    9: ("-36188894799/29900288", ["1", "0", "74376/5309", "5258730/58399", "25667658/58399", "0", "0", "0",
                                  "817181091/116798", "134340290/58399"], [3, 4, 5, 8, 10, 24],
        ["1248777/467192", "5312223/3737536", "3466395/3737536", "335192013/7475072", "2086155/934384",
         "23571/1868768"], "-5515759/22425216"),
}
for n, (const, Avec, active, mu, qlow) in DATA.items():
    A, cons, obj = enum_lp(n)
    const = Fr(const); mu = [Fr(t) for t in mu]; x = [Fr(t) for t in Avec[2:-2]]
    recon = [sum(mu[i] * cons[a][1][c] for i, a in enumerate(active)) + (const if c == 0 else 0)
             for c in range(len(obj))]
    dual = recon == obj and all(t >= 0 for t in mu)
    ev = lambda aff: aff[0] + sum(aff[i + 1] * x[i] for i in range(len(x)))
    primal = all(ev(c[1]) >= 0 for c in cons) and ev(obj) == const
    full = [ev(A[j]) for j in range(n + 1)] == [Fr(t) for t in Avec]
    lower = Fr(4, 3 ** n) * const
    good = dual and primal and full and lower == Fr(qlow)
    ok &= good
    print(f'n={n}: dual identity {dual}, primal feasible {primal}, primal vector consistent {full}, '
          f'bound 4/3^n * const = {lower} ({float(lower):+.5f}) matches {lower == Fr(qlow)}')
print('PASS' if ok else 'FAIL')
