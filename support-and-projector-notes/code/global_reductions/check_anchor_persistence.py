"""Doubly balanced encoding keeps anchor obstructions: at r = 1, the column flag word y = 1...1 spans an invariant
block of D (no coupling to any other input column flag word), only branch 0 reaches (row flags 0...0, column flags y),
and the diagonal weights are 1 (branch 0) and 2^{-3r} (other branches). Exhaustive over all 3^8 column flag words."""
import itertools
LAB = [(a, b) for a in range(3) for b in range(3)]
ELL = [l for l in LAB if l != (0, 0)]
r = 1
flags = [tuple(((l[0] * lam[0] + l[1] * lam[1]) % 3) for l in ELL for _ in range(r)) for lam in LAB]
y = (1,) * (8 * r); c0 = (0,) * (8 * r)

def m(a, c, b, d):   # <a, c| M |b, d>, rows a, b; columns c, d
    v = 1.0
    for s in range(len(a)):
        v *= (a[s] == b[s] and c[s] == d[s]) - 0.5 * (a[s] == c[s] and b[s] == d[s])
        if v == 0:
            return 0.0
    return v

words = list(itertools.product(range(3), repeat=8 * r))
Doff = sum(1 for a in range(9) for b in range(9) for z in words
           if m(flags[a], y, flags[b], z) != 0 and (a != b or z != y))
Goff = sum(1 for a in range(9) for z in words if m(c0, y, flags[a], z) != 0 and (a != 0 or z != y))
diag = [m(flags[a], y, flags[a], y) for a in range(9)]
print('D couplings off the block:', Doff, ' MV entries off branch 0:', Goff, ' diagonal weights:', diag)
print('PASS' if Doff == 0 and Goff == 0 and diag == [1.0] + [2.0 ** (-3 * r)] * 8 else 'FAIL')
