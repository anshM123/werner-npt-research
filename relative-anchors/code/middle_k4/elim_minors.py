"""Leading-principal-minor certificate for one gauge-fixed torus block (exact, python-flint).

Input: blocks_sym.pkl (build_sym.py).  The gauge-fixed block Mk (entries in Q(s1, s2, s3), s0 = 1) is reordered by the
chosen order and eliminated symmetrically over the rational function field.  The r-th pivot is D_r/D_{r-1}, D_r the
leading r x r minor; D_r = D_{r-1} * pivot_r is kept as a reduced fraction N_r/Q_r.  The certificate for step r is:
Q_r and N_r each have coefficients of one sign, the same sign (so D_r > 0 on the open orthant).  If every step passes,
all leading principal minors are positive at every point of the open orthant, so Mk > 0 there (Sylvester).  Entries
are homogeneous in (s0, ..., s3), so coefficient signs with s0 = 1 are those of the homogeneous polynomials.
Usage: python elim_minors.py j T order [maxsteps]
   T: '-' for the letter-1 block, else e.g. 0 or 01 or 012;  order: 'nat', 'rev', 'vacfirst', 'vaclast', 'col', or a
   comma-separated permutation of 0..n-1 (positions in the gauge-fixed list)."""
import sys, time, pickle, os
import flint
from flintfield import Field, F as FF

j = int(sys.argv[1])
T = () if sys.argv[2] == '-' else tuple(int(c) for c in sys.argv[2])
ordspec = sys.argv[3] if len(sys.argv) > 3 else 'nat'
maxsteps = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
D = pickle.load(open('blocks_sym.pkl', 'rb'))
K = Field(D['names'])
blk = D['blocks'][(j, T)]
S = blk['states']
keep = blk['keep']
n = len(keep)
ctx = K.ctx


def fe(nd):
    return FF(K, ctx.from_dict(nd[0]), ctx.from_dict(nd[1]))


M = [[fe(blk['M'][a][b]) for b in keep] for a in keep]
labs = [blk['labels'][a] for a in keep]
vac = (0, 0, 0, 0)
isvac = [S[a][0] == vac for a in keep]
if ordspec == 'nat':
    order = list(range(n))
elif ordspec == 'rev':
    order = list(range(n))[::-1]
elif ordspec == 'vacfirst':
    order = [i for i in range(n) if isvac[i]] + [i for i in range(n) if not isvac[i]]
elif ordspec == 'vaclast':
    order = [i for i in range(n) if not isvac[i]] + [i for i in range(n) if isvac[i]]
elif ordspec == 'col':
    order = sorted(range(n), key=lambda i: (isvac[i], S[keep[i]][1], S[keep[i]][0]))
else:
    order = [int(x) for x in ordspec.split(',')]
assert sorted(order) == list(range(n))
A = [[M[a][b] for b in order] for a in order]
print(f'j = {j}, T = {T}: gauge-fixed size {n}, order {ordspec}: {[labs[i] for i in order]}', flush=True)


def sgn(p):
    cs = p.coeffs()
    if all(c > 0 for c in cs):
        return 1
    if all(c < 0 for c in cs):
        return -1
    return 0


t0 = time.time()
Dr = K(1)
ok = True
res = []
for p in range(min(n, maxsteps)):
    pv = A[p][p]
    Dr = Dr * pv
    sn, sd = sgn(Dr.n), sgn(Dr.d)
    good = sn != 0 and sn == sd
    nneg = sum(1 for c in Dr.n.coeffs() if c < 0)
    psn, psd = sgn(pv.n), sgn(pv.d)
    res.append((p + 1, good, Dr.n.total_degree(), len(Dr.n), nneg))
    print(f'  D_{p + 1}: num deg {Dr.n.total_degree()}, {len(Dr.n)} terms, {nneg} negative; den deg {Dr.d.total_degree()}; '
          f'minor certified > 0: {good}; pivot itself one-signed: {psn != 0 and psn == psd}  [{time.time() - t0:.1f}s]',
          flush=True)
    ok &= good
    if pv == 0:
        break
    for i in range(p + 1, n):
        if A[i][p] == 0:
            continue
        f = A[i][p] / pv
        for l in range(i, n):
            if A[p][l] != 0:
                A[i][l] = A[i][l] - f * A[p][l]
                if l != i:
                    A[l][i] = A[i][l]
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok and len(res) == n else 'FAIL')
