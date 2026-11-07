"""Exact (sympy) checks of the scalar facts used in NOTES.md, Sections 2-4."""
import sympy as sp
t, s, x = sp.symbols('t s x', positive=True)
K = t + sp.Rational(1, 2)
D = 2 * K + 1
ok = True
def check(name, cond):
    global ok
    ok &= bool(cond)
    print(('ok   ' if cond else 'FAIL ') + name)
# Section 3: antisymmetric hard bound and the d-channel loss
f = s ** 2 / (4 * (2 + 3 * s / 2))
check('s^2/(4(2 + 3s/2)) is increasing in s > 0', sp.simplify(sp.diff(f, s) - s * (3 * s + 8) / (2 * (3 * s + 4) ** 2)) == 0)
check('its value at s = 2/(3K+1) is 1/((3K+1)(6K+5))', sp.simplify(f.subs(s, 2 / (3 * K + 1)) - 1 / ((3 * K + 1) * (6 * K + 5))) == 0)
h = sp.Rational(1, 2) - 1 / ((3 * K + 1) * (6 * K + 5))
delta = t / (8 * K ** 2 * (K + 1) ** 3)
N = sp.expand(sp.simplify((h - delta - sp.Rational(4, 9)) * 144 * K ** 2 * (K + 1) ** 3 * (3 * K + 1) * (6 * K + 5)))
claimed = (8 * (t - 1) ** 2 + sp.Rational(11, 2) + t / 4 + 1171 * t ** 2 + 4099 * t ** 3 + 5296 * t ** 4 + 3388 * t ** 5
           + 1104 * t ** 6 + 144 * t ** 7)
check('(h - delta - 4/9) 144 K^2 (K+1)^3 (3K+1)(6K+5) = 8(t-1)^2 + 11/2 + t/4 + 1171 t^2 + ... > 0', sp.expand(N - claimed) == 0)
check('x^2/(K+x) <= 1/(K+1) on [0,1]: (K + x) - (K+1) x^2 = K(1 - x^2) + x(1 - x) >= 0',
      sp.expand((K + x) - (K + 1) * x ** 2 - (K * (1 - x ** 2) + x * (1 - x))) == 0)
check('Popoviciu: Var of f on [0,1] <= ((1/K - 1/(K+1))/2)^2 = 1/(4 K^2 (K+1)^2)',
      sp.simplify(((1 / K - 1 / (K + 1)) / 2) ** 2 - 1 / (4 * K ** 2 * (K + 1) ** 2)) == 0)
# Section 2: the constant R* and the soft lower bound
tau = t / (t + 1)
check('t/K >= tau: t/K - tau = t/(2 K (t+1)) > 0', sp.simplify(t / K - tau - t / (2 * K * (t + 1))) == 0)
check('t(t^2 + t + 1/2)/(2 K^2 (t+1)) - tau/2 = t/(8 K^2 (t+1)) > 0',
      sp.simplify(t * (t ** 2 + t + sp.Rational(1, 2)) / (2 * K ** 2 * (t + 1)) - tau / 2 - t / (8 * K ** 2 * (t + 1))) == 0)
check('4 K^4 - t D (t D + 1) = 1/4', sp.expand(4 * K ** 4 - t * D * (t * D + 1)) == sp.Rational(1, 4))
# Section 4: shells
check('1 - sqrt(2/3) > 1/6', sp.N(1 - sp.sqrt(sp.Rational(2, 3))) > sp.Rational(1, 6))
check('t/(2(K+1)) >= tau/3: difference t^2/(3(t+1)(2t+3)) >= 0',
      sp.simplify(t / (2 * (K + 1)) - tau / 3 - t ** 2 / (3 * (t + 1) * (2 * t + 3))) == 0)
check('min(tau/12, 2 tau/9) = tau/12', sp.simplify(sp.Min(tau / 12, 2 * tau / 9) - tau / 12) == 0)
print('PASS' if ok else 'FAIL')
