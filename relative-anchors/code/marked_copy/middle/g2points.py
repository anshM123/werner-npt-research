"""An independent check of the interpolation (not needed by the proof): at the four test points of stage B, all
outside the grid, det Q * Den_c computed directly from the form (every variable specialised: a 10 x 10 rational
determinant) must equal the value of N~ accumulated in stage B (OUT/testpoints.pkl).
python g2points.py OUT"""
import sys, os, pickle
import flint
import g2core, g2slices

out = sys.argv[1]
S = g2core.Slicer(check=False, log=lambda m: None)
Dc = g2slices.den_c(S)
with open(os.path.join(out, 'testpoints.pkl'), 'rb') as fh:
    tp = pickle.load(fh)
ok = True
for pt, (vp, vq) in tp:
    x, y, z, s = (flint.fmpq(p, q) for p, q in pt)
    vals = [x + 4, y + 4, z + 7, s]
    A = [[(S.M[a, b].n(*vals) / S.M[a, b].d(*vals)) if S.M[a, b] != 0 else flint.fmpq(0) for b in S.keep]
         for a in S.keep]
    det = flint.fmpq_mat(len(A), len(A), [v for row in A for v in row]).det()
    lhs = det * Dc(*vals)
    rhs = flint.fmpq(vp, vq)
    good = lhs == rhs
    ok &= good
    print(f'(x, y, z, s) = {tuple(str(v) for v in (x, y, z, s))}: det Q * Den_c = {float(lhs):.6e}, N~ = '
          f'{float(rhs):.6e}, equal: {good}')
print('PASS' if ok else 'FAIL')
