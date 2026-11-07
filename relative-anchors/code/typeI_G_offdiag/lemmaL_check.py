"""Symbolic checks for Lemma L (local multiplier):  L(s) = (K+s)^2/((ka+s)(kb+s)) - sigma/(k0+s),
K = sigma + 1/2, ka = K - a/2, kb = K - b/2, k0 = K - (a+b)/2, a, b >= 0, a + b <= sigma, s >= 0.
With d = K + s, A = a/2, B = b/2, S' = A + B, P = A B:
 (i)  (1 - L) (d-A)(d-B)(d-S') = d (sigma - S')(d - S') + P (sigma + d - S')        [identity]
 (ii) 1/(1-q/2)^2 - 1/(1-q) = -q^2/(4 (1-q)(1-q/2)^2)  and  q^2/(4(1-q)(1-q/2)^2) at q = 1/2 equals 2/9, and is
      increasing on [0, 1/2] (derivative numerator has positive coefficients after q -> q).
 (iii) L = 1 + (a/2) g_a + (b/2) g_b + (ab/4) g_a g_b - sigma g_0   [partial fractions used for E[L Z°]]."""
import sympy as sp
sig, a, b, s, q = sp.symbols('sigma a b s q', nonnegative=True)
K = sig + sp.Rational(1, 2); ka, kb, k0 = K - a/2, K - b/2, K - (a+b)/2
L = (K+s)**2/((ka+s)*(kb+s)) - sig/(k0+s)
d = K + s; A_, B_ = a/2, b/2; Sp = A_ + B_; P = A_*B_
lhs = (1 - L)*(d-A_)*(d-B_)*(d-Sp)
rhs = d*(sig - Sp)*(d - Sp) + P*(sig + d - Sp)
print('(i) identity:', sp.simplify(lhs - rhs) == 0)
e1 = 1/(1-q/2)**2 - 1/(1-q) + q**2/(4*(1-q)*(1-q/2)**2)
print('(ii) identity:', sp.simplify(e1) == 0, '; value at q=1/2:', sp.nsimplify((q**2/(4*(1-q)*(1-q/2)**2)).subs(q, sp.Rational(1, 2))))
h = q**2/(4*(1-q)*(1-q/2)**2)
num = sp.factor(sp.numer(sp.together(sp.diff(h, q))))
print('     derivative numerator:', num, '(positive on (0, 1/2]: q(... ) with (1-q)>0 factors)')
ga, gb, g0 = 1/(ka+s), 1/(kb+s), 1/(k0+s)
print('(iii) partial fractions:', sp.simplify(L - (1 + a/2*ga + b/2*gb + a*b/4*ga*gb - sig*g0)) == 0)
