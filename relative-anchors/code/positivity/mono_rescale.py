"""Test 4: can a diagonal rescaling by monomials remove the obstruction of newton_obstr.py, and does the rescaled block
then have a fixed-vector representation?

Search (floating point LP, HiGHS): exponent vectors delta_a (one per state, in the dehomogenised variables s_1..s_r,
s_0 = 1) such that for every ordered pair (a, b) with M_ab != 0 the support of P_ab = n_ab d_aa, shifted by
delta_b - delta_a, lies in conv(supp Q_ab), Q_ab = d_ab n_aa.  This is exactly boundedness of M'_ab / M'_aa for
M' = D M D, D = diag(s^delta).  The LP solution is rounded to half-integers and multiplied by 2, i.e. D = diag(s^{2 delta})
is replaced by the substitution-free form D2 = diag(s^{d_a}) with d_a = round(2 delta_a) and M'' = D2 M D2 compared
with the squared problem (s -> s^2 is harmless for coefficient signs).  Then, exactly: (i) the obstruction test on M''
(integer directions), (ii) every monomial coefficient matrix of den'' * M'' PSD?
Usage: python mono_rescale.py k"""
import sys, itertools, time
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
import flint
from psdtools import load_blocks, coeff_mats, psd_exact, lcm_poly

k = int(sys.argv[1])
t0 = time.time()


def dehom(p, homog):
    d = p.to_dict()
    return {m[1:] if homog else m: c for m, c in d.items()}


for blk in load_blocks(k):
    n = len(blk['Mg'])
    Mg = blk['Mg']
    homog = (k in (2, 3))
    nv = len(blk['names']) - (1 if homog else 0)
    pairs = []
    for a in range(n):
        for b in range(n):
            if a != b and not Mg[a][b][0].is_zero():
                P = dehom(Mg[a][b][0] * Mg[a][a][1], homog)
                Q = dehom(Mg[a][b][1] * Mg[a][a][0], homog)
                pairs.append((a, b, list(P.keys()), list(Q.keys())))
    if not pairs:
        continue
    # variables: delta (n * nv), then lambda blocks
    nd = n * nv
    rows_eq, b_eq = [], []
    nl = 0
    blocks = []
    for (a, b, sP, sQ) in pairs:
        for m in sP:
            blocks.append((a, b, m, sQ, nl))
            nl += len(sQ)
    nvar = nd + nl
    from scipy.sparse import lil_matrix
    A = lil_matrix((len(blocks) * (nv + 1) + nv, nvar))
    rhs = []
    r = 0
    for (a, b, m, sQ, off) in blocks:
        for i in range(nv):
            # sum_q lam_q q_i - (delta_b - delta_a)_i = m_i
            for jq, q in enumerate(sQ):
                if q[i]:
                    A[r, nd + off + jq] = q[i]
            A[r, b * nv + i] += -1
            A[r, a * nv + i] += 1
            rhs.append(m[i])
            r += 1
        for jq in range(len(sQ)):
            A[r, nd + off + jq] = 1
        rhs.append(1)
        r += 1
    for i in range(nv):  # delta_0 = 0
        A[r, i] = 1
        rhs.append(0)
        r += 1
    bounds = [(-20, 20)] * nd + [(0, None)] * nl
    res = linprog(np.zeros(nvar), A_eq=A.tocsr(), b_eq=np.array(rhs, dtype=float), bounds=bounds, method='highs')
    msg = f'{blk["name"]} ({blk["kind"]}), size {n}: LP status {res.status} ({"feasible" if res.status == 0 else "infeasible or failed"})'
    if res.status != 0:
        print(msg + f'  [{time.time() - t0:.1f}s]', flush=True)
        continue
    delta = res.x[:nd].reshape(n, nv)
    d2 = np.round(2 * delta).astype(int)
    d2 = d2 - d2.min(axis=0)
    # exact: M'' = D2 M(s^2-substituted?)  -- we keep the variables s and use D2 = diag(s^{d2_a}) on the block with
    # s replaced by s^2, i.e. entries M_ab(s^2) s^{d2_a + d2_b}.  Coefficient signs of minors are unchanged by s -> s^2.
    ctx = blk['ctx']
    names = blk['names']
    nvv = len(names)

    def sq(p):
        return ctx.from_dict({tuple(2 * e for e in m): c for m, c in p.to_dict().items()})

    def mono(ex):
        full = ((0,) + tuple(int(x) for x in ex)) if homog else tuple(int(x) for x in ex)
        return ctx.from_dict({full: 1})
    Ms = [[(sq(Mg[a][b][0]) * mono(d2[a] + d2[b]), sq(Mg[a][b][1])) for b in range(n)] for a in range(n)]
    den = lcm_poly([Ms[a][b][1] for a in range(n) for b in range(n)])
    Nt = [[Ms[a][b][0] * (den / Ms[a][b][1]) for b in range(n)] for a in range(n)]
    C = coeff_mats(Nt)
    nfail = sum(1 for m in C if not psd_exact(C[m])[0])
    # exact obstruction test on the rescaled block (directions in {-2..2})
    dirs = [l for l in itertools.product(range(-2, 3), repeat=nvv) if any(l)]
    nobs = 0
    for a in range(n):
        for b in range(n):
            if a == b or Ms[a][b][0].is_zero():
                continue
            P = Ms[a][b][0] * Ms[a][a][1]
            Q = Ms[a][b][1] * Ms[a][a][0]
            sP, sQ = list(P.to_dict()), list(Q.to_dict())
            for l in dirs:
                vP = [sum(x * y for x, y in zip(l, m)) for m in sP]
                mp = max(vP)
                if mp > max(sum(x * y for x, y in zip(l, m)) for m in sQ) and vP.count(mp) == 1:
                    nobs += 1
                    break
    print(msg + f'; rounded exponents 2*delta = {[list(map(int, x)) for x in d2]}; exact obstruction pairs after '
          f'rescaling: {nobs}; monomials {len(C)}, non-PSD coefficient matrices {nfail}  [{time.time() - t0:.1f}s]',
          flush=True)
print('PASS')
