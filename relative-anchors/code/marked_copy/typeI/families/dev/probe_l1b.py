"""Development probe: the rows without the T' letter against the letter-1 form, with B scaled by b (b = 0: the A-parts
only).  Prints the largest |difference| / |entry| ratio pattern."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import tlib, mlib

def run(q0, t, i, tp, sg, b):
    sg, b = Fr(sg), Fr(b)
    k = tp + q0 + t + 2 * i + 1
    j, tau = q0 + tp + i, tp
    n, _, Q = tlib.reduced(tlib.Calc(k, j, tau, i, 'MU', 'sym', sigma=sg, bscale=b))
    n1, _, Q1 = mlib.reduced(mlib.Calc(k - tp, j - tp, i, sigma=sg, bscale=b))
    ix = {nm: a for a, nm in enumerate(n)}
    ix1 = {nm: a for a, nm in enumerate(n1)}
    rat = {}
    for a in n1:
        for c in n1:
            x = Q[ix[a]][ix[c]]; y = Q1[ix1[a]][ix1[c]]
            x = Fr(int(x.p), int(x.q)); y = Fr(int(y.p), int(y.q))
            if y != 0:
                rat[(a, c)] = x / y
            elif x != 0:
                rat[(a, c)] = 'inf'
    vals = sorted(set(str(v) for v in rat.values()))
    print(f'b={b}: distinct ratios tlib/mlib on the letter-1 vectors: {vals[:12]}{" ..." if len(vals) > 12 else ""}')
    for key in [('O0', 'O0'), ('A00', 'A00'), ('B000', 'B000'), ('G0', 'G0'), ('O0', 'G1'), ('A10', 'B101')]:
        print('   ', key, rat.get(key))

for b in ('0', '1'):
    run(3, 3, 2, 2, '1/3', b)
