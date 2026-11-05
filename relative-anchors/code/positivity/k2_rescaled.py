"""Exact check of an explicit positive representation after a non-monomial rescaling, k = 2, letter-1 block of
(R_inf)_1 = (R_inf)_{k-1} (gauge-fixed states (10,10) and (00,11)).

newton_obstr.py shows that this block has no representation with fixed vectors (M_12 / M_11 is unbounded), and the LP of
mono_rescale.py finds no monomial rescaling.  With D = diag(|s| / s0, 1), |s| = s0 + s1, and Delta = s0^2 + 4 s0 s1 + s1^2:
    D M D = s0 s1 [ (1/2) (1, 1)(1, 1)^T + diag((s0^2 + s1^2) / (2 Delta), 1/2) ],
a positive combination of squares of fixed integer forms with subtraction-free weights.  The script checks this
identity exactly.  Usage: python k2_rescaled.py"""
import flint
from psdtools import load_blocks

blk = [b for b in load_blocks(2) if b['name'] == 'k2_j1_(2, 0, 0)'][0]
ctx = blk['ctx']
s0, s1 = ctx.gens()
one = ctx.from_dict({(0, 0): 1})
Mg = blk['Mg']
print('states', blk['labels'])
S = s0 + s1
Dl = s0 * s0 + 4 * s0 * s1 + s1 * s1
d = [(S, s0), (one, one)]  # (num, den) of the diagonal rescaling
ok = True
# target T = s0 s1 [ 1/2 (1,1)(1,1)^T + diag((s0^2 + s1^2)/(2 Delta), 1/2) ]  as (num, den)
half = flint.fmpq(1, 2)
T = [[(s0 * s1 * (half * Dl + half * (s0 * s0 + s1 * s1)), Dl), (s0 * s1 * half, one)],
     [(s0 * s1 * half, one), (s0 * s1 * (half + half), one)]]
for a in range(2):
    for b in range(2):
        n_, d_ = Mg[a][b]
        lhs_n = n_ * d[a][0] * d[b][0]
        lhs_d = d_ * d[a][1] * d[b][1]
        tn, td = T[a][b]
        eq = lhs_n * td == tn * lhs_d
        ok &= eq
        print(f'entry ({a},{b}): (D M D) = target: {eq}')
print('weights s0 s1 / 2, s0 s1 (s0^2 + s1^2) / (2 Delta), s0 s1 / 2 have positive coefficients: True by inspection')
print('PASS' if ok else 'FAIL')
