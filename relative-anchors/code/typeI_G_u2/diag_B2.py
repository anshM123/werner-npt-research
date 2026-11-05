import pickle, flint, itertools
from fractions import Fraction as Fr
import build_B as bb
num, den = bb.load()
th, mu, la = bb.th, bb.mu, bb.la
H = flint.fmpq(2, 5)
sig = 1 - mu; st = th * sig; sU = [mu * la, mu * (1 - la)]; h = H * (1 - th); one = bb.one
Oinv = [[one, 0 * one, 0 * one], [0 * one, st * (sU[0] - sU[0] ** 2), -st * sU[0] * sU[1]], [0 * one, -st * sU[0] * sU[1], st * (sU[1] - sU[1] ** 2)]]
Nt = [[num[i][j] - h * den * Oinv[i][j] for j in range(3)] for i in range(3)]
cands = {'z3': Nt[1][1], 'z4': Nt[2][2], 'z3z4': Nt[1][1] * Nt[2][2] - Nt[1][2] ** 2, 'rz3': Nt[0][0] * Nt[1][1] - Nt[0][1] ** 2,
         'rz4': Nt[0][0] * Nt[2][2] - Nt[0][2] ** 2}
pts = []
g = [Fr(i, 6) for i in range(7)]
for face in range(3):
    for val in (Fr(0), Fr(1)):
        for a in g:
            for b in g:
                x = [a, b]; x.insert(face, val); pts.append(tuple(x))
def ev(p, x):
    return sum(Fr(int(c.p), int(c.q)) * x[0] ** int(m[0]) * x[1] ** int(m[1]) * x[2] ** int(m[2]) for m, c in p.to_dict().items())
for nm, p in cands.items():
    q, rem = bb.strip(p)
    zeros = sorted(set(x for x in pts if ev(q, x) <= 0))
    print(nm, 'degrees', q.degrees(), 'removed', rem, 'nonpositive boundary samples:', [tuple(str(v) for v in x) for x in zeros][:12], len(zeros))
