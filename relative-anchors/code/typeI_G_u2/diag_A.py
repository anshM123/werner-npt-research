import flint, itertools, time
from fractions import Fraction as Fr
import build_A as ba
FA, den = ba.build(Fr(2, 5))
q, rem = ba.strip(FA)
F = flint.fmpq
g = [F(i, 8) for i in range(9)]
neg = []; zero = []
mn = None
t0 = time.time()
for th_ in g:
    for ph_ in g:
        if th_ + ph_ > 1:
            continue
        for mu_ in g:
            for la_ in g:
                on_bd = th_ in (0,) or ph_ in (0,) or th_ + ph_ == 1 or mu_ in (0, 1) or la_ in (0, 1)
                v = q(th_, ph_, mu_, la_)
                if v < 0: neg.append((th_, ph_, mu_, la_))
                if v == 0: zero.append((th_, ph_, mu_, la_))
                # normalised value
                if mn is None or v < mn[0]:
                    mn = (v, (th_, ph_, mu_, la_))
print('negatives', [tuple(str(x) for x in p) for p in neg][:20], len(neg))
print('zeros', [tuple(str(x) for x in p) for p in zero][:40], len(zero))
print('min', float(mn[0]), [str(x) for x in mn[1]], f'{time.time() - t0:.1f}s')
