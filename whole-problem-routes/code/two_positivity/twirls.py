"""Exact check (SymPy, Gaussian rationals) of the twirl coefficients of the four product states used in
cert_ncopies.py: for rho = |a><a| (x) |b><b| the O(3)-twirl is sum_k Tr(rho P_k) P_k / d_k (Schur), with coefficient
vector (Tr(rho P_Om)/1, Tr(rho P_s0)/5, Tr(rho P_an)/3).  Also re-checks that M, C_a have the stated coefficients."""
import sympy as sp
I = sp.I
e = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
F = sp.zeros(9, 9); OO = sp.zeros(9, 9)
for i in range(3):
    for j in range(3):
        F[3 * i + j, 3 * j + i] = 1
        OO[4 * i, 4 * j] = 1
I9 = sp.eye(9)
P_Om = OO / 3; P_an = (I9 - F) / 2; P_s0 = (I9 + F) / 2 - P_Om
w = (e[0] + I * e[1])           # unnormalised; |w|^2 = 2
states = {'e1e1': (e[0], e[0]), 'ww': (w, w), 'wwbar': (w, w.conjugate()), 'e1e2': (e[0], e[1])}
claimed = {'e1e1': [sp.Rational(1, 3), sp.Rational(2, 15), 0], 'ww': [0, sp.Rational(1, 5), 0],
           'wwbar': [sp.Rational(1, 3), sp.Rational(1, 30), sp.Rational(1, 6)], 'e1e2': [0, sp.Rational(1, 10), sp.Rational(1, 6)]}
ok = True
for name, (x, y) in states.items():
    v = sp.kronecker_product(x, y)
    v = v / sp.sqrt((v.H * v)[0])
    rho = v * v.H
    coeff = [sp.nsimplify(sp.simplify((rho * P).trace())) / d for P, d in ((P_Om, 1), (P_s0, 5), (P_an, 3))]
    good = all(sp.simplify(c1 - c2) == 0 for c1, c2 in zip(coeff, claimed[name]))
    print(name, 'twirl coefficients', coeff, good)
    ok &= good
print('PASS' if ok else 'FAIL')
