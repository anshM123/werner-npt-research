"""Rigorous interval-arithmetic verification of the two reduced inequalities behind (C0) (bottom condition (R_inf)_0,
letter-1 block), every k, every weight vector (relative-anchors/README.md, Theorem B0).

Regime A (x1 + x2 <= 3/4, x2 <= x1):     fA = l1 + 2 + x1 l2 (l2 + 2)  <=  (2 + gA)^2,
   l_i = 1/(1 - x_i), gA = (l1 - 1)/2 + (l2 - 1)/2 + (1 - x1 - x2)/2.
Regime B (x1 + x2 > 3/4; s = x1 + x2 in (3/4, 1], b = x2 in (0, s/2], a = s - b, eps = 1 - s):
   U2 = max over the vertices (q1, q2) of {0 <= q1 <= 1-a, 0 <= q2 <= 1-b, q1 + q2 <= 1} of
        c1 q1 + c2 q2 + ct (1 - q1 - q2), plus (1/2) l1 l2,
   c1 = l1^2 + l2^2,  c2 = l2^2 + l2 min(l1, 1/(4b)) + (4/3) l2^2 eps,  ct = lt^2 + lt l1 + (4/3) lt^2 (b + eps),
   lt = 1/(1 - eps);  need U2 <= (2 + gB)^2 with gB = (l1 - 1)/2 + (l2 - 1)/2 + eps/2.
   The corner a >= 0.99 is excluded here and handled analytically (ratio <= 0.1, see the notes).
Each box of an adaptive subdivision is evaluated in interval arithmetic (mpmath.iv); a box passes if the upper bound of
LHS - RHS is < 0.  Prints PASS if every box passes."""
from mpmath import iv, mpf
iv.dps = 30
def imin(x, y): return iv.mpf([min(x.a, y.a), min(x.b, y.b)])
def imax(x, y): return iv.mpf([max(x.a, y.a), max(x.b, y.b)])

def checkA(x1, x2):
    l1 = 1 / (1 - x1); l2 = 1 / (1 - x2)
    f = l1 + 2 + x1 * l2 * (l2 + 2)
    g = (l1 - 1) / 2 + (l2 - 1) / 2 + (1 - x1 - x2) / 2
    return (f - (2 + g) ** 2).b < 0

def checkB(s, b):
    a = s - b; eps = 1 - s
    l1 = 1 / (1 - a); l2 = 1 / (1 - b); lt = 1 / (1 - eps)
    c1 = l1 ** 2 + l2 ** 2
    m = imin(l1, 1 / (4 * b)) if b.a > 0 else l1     # min(l1, 1/(4b)) <= l1 (we only need an upper bound)
    c2 = l2 ** 2 + l2 * m + iv.mpf(4) / 3 * l2 ** 2 * eps
    ct = lt ** 2 + lt * l1 + iv.mpf(4) / 3 * lt ** 2 * (b + eps)
    verts = [(iv.mpf(0), iv.mpf(0)), (1 - a, iv.mpf(0)), (iv.mpf(0), 1 - b), (1 - a, a), (b, 1 - b)]
    best = None
    for q1, q2 in verts:
        v = c1 * q1 + c2 * q2 + ct * (1 - q1 - q2)
        best = v if best is None else imax(best, v)
    g = (l1 - 1) / 2 + (l2 - 1) / 2 + eps / 2
    return (best + l1 * l2 / 2 - (2 + g) ** 2).b < 0

def run(check, box, inside, maxdepth=22):
    """box = ((lo1, hi1), (lo2, hi2)); inside(lo1,hi1,lo2,hi2) -> 'out' | 'in' | 'partial'. Subdivide until pass."""
    stack = [(box, 0)]; n_ok = 0; worst_depth = 0
    while stack:
        ((u0, u1), (v0, v1)), d = stack.pop()
        st = inside(u0, u1, v0, v1)
        if st == 'out': continue
        X = iv.mpf([u0, u1]); Yv = iv.mpf([v0, v1])
        if check(X, Yv):   # boxes meeting the region are evaluated whole (formulas are defined on the whole box)
            n_ok += 1; worst_depth = max(worst_depth, d); continue
        if d >= maxdepth:
            return False, (u0, u1, v0, v1)
        um = (u0 + u1) / 2; vm = (v0 + v1) / 2
        stack += [(((u0, um), (v0, vm)), d + 1), (((um, u1), (v0, vm)), d + 1),
                  (((u0, um), (vm, v1)), d + 1), (((um, u1), (vm, v1)), d + 1)]
    return True, (n_ok, worst_depth)

def insideA(u0, u1, v0, v1):
    # region: 0 < x2 <= x1, x1 + x2 <= 3/4 ; boxes in (x1, x2)
    if v0 > u1 or u0 + v0 > mpf(3) / 4 or v1 <= 0: return 'out'
    if v1 <= u0 and u1 + v1 <= mpf(3) / 4 and v0 > 0: return 'in'
    return 'partial'

def insideB(u0, u1, v0, v1):
    # (s, b): 3/4 < s <= 1, 0 < b <= s/2, a = s - b <= 0.99
    if u1 <= mpf(3) / 4 or v1 <= 0 or v0 > u1 / 2 or u0 - v1 > mpf('0.99'): return 'out'
    if u0 > mpf(3) / 4 and v0 > 0 and v1 <= u0 / 2 and u1 - v0 <= mpf('0.99'): return 'in'
    return 'partial'

if __name__ == '__main__':
    import time
    t0 = time.time()
    okA, infoA = run(checkA, ((mpf('1e-9'), mpf(3) / 4), (mpf('1e-9'), mpf(3) / 8)), insideA, maxdepth=14)
    print('Regime A:', 'PASS' if okA else 'FAIL', infoA, f'[{time.time()-t0:.0f}s]', flush=True)
    okB, infoB = run(checkB, ((mpf(3) / 4, mpf(1)), (mpf('1e-9'), mpf(1) / 2)), insideB, maxdepth=14)
    print('Regime B (a <= 0.99):', 'PASS' if okB else 'FAIL', infoB, f'[{time.time()-t0:.0f}s]', flush=True)
