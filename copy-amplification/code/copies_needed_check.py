"""Exact check of n_min(alpha) = m (kappa_PPT - kappa_0)/(kappa_PPT - kappa(alpha)) against its closed forms."""
from fractions import Fraction as Fr


def nmin(alpha, d, m, a0):
    kp = Fr(d - 1, d + 1); k0 = (1 - a0) / (1 + a0); k = (1 - alpha) / (1 + alpha)
    return m * (kp - k0) / (kp - k)


for a in (Fr(35, 100), Fr(34, 100), Fr(5, 11), Fr(2, 5)):
    assert nmin(a, 3, 3, Fr(1, 2)) == (1 + a) / (3 * a - 1)
    print('d=3, m=3, alpha =', a, ': n_min =', nmin(a, 3, 3, Fr(1, 2)))
for d in (3, 4, 5, 10, 50):
    for a in (Fr(1, d) + Fr(1, 60), Fr(1, d) + Fr(1, 7)):
        assert nmin(a, d, 2, Fr(1, 2)) == 2 * (d - 2) * (1 + a) / (3 * (a * d - 1))
print('closed forms confirmed')
