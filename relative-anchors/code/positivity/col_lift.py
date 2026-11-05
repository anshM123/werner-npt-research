"""Column groups: the states of a gauge-fixed block with one fixed column string y.  On a column group the partial
transpose does nothing, so the B-part is an ordinary Gram matrix:
    Y'[(a,y),(b,y)] = A'[(a,y),(b,y)] - <MV'(b,y), D'^{-1} MV'(a,y)>,
and the column group is the Schur complement onto x of the lift  L = [[A'_col, W^T], [W, D']], W[z, a] = MV'(a, y)[z].
Checks: (1) exact at a rational point: Schur complement of L = column group, and L PSD; (2) symbolic: obstruction,
Farkas against monomial rescalings, LP rescaling + monomial-PSD test.
Usage: python col_lift.py k"""
import sys, time
import flint
from ylib import Builder, Ring, point_ring, wcol
from flintfield import Field
from psdtools import load_blocks
from analyze import obstruction, mono_psd, farkas, rescale_lp, apply_rescale
from lift_test_lib import Dprime_blocks, psd_point

k = int(sys.argv[1])
t0 = time.time()
summary = []
for blk in load_blocks(k):
    states = blk['states']
    groups = {}
    for st in states:
        if any(c != 0 for c in st[0]):
            groups.setdefault(st[1], []).append(st)
    for y, G in sorted(groups.items()):
        if len(G) < 2:
            continue

        def build(Bd, zero):
            Dp = Dprime_blocks(Bd, k)
            idx = {('x', st): i for i, st in enumerate(G)}
            vecs = {st: Bd.MV(st[0], st[1]) for st in G}
            keys = {wcol(z) for st in G for z in vecs[st]}
            for key in sorted(keys):
                for z in Dp[key][0]:
                    idx[('z', z)] = len(idx)
            n = len(idx)
            M = [[zero] * n for _ in range(n)]
            for X1 in G:
                for X2 in G:
                    M[idx[('x', X1)]][idx[('x', X2)]] = Bd.A(X1[0], X1[1], X2[0], X2[1])
            for key in keys:
                ys, Bm = Dp[key]
                for p, z1 in enumerate(ys):
                    for q, z2 in enumerate(ys):
                        M[idx[('z', z1)]][idx[('z', z2)]] = Bm[p][q]
            for st in G:
                for z, c in vecs[st].items():
                    M[idx[('z', z)]][idx[('x', st)]] = c
                    M[idx[('x', st)]][idx[('z', z)]] = c
            return M, len(G)
        w = [flint.fmpq(p, q) for p, q in [(2, 3), (5, 4), (3, 7), (7, 5)][:k]]
        sp = [x * x for x in w]
        Bp = Builder(point_ring(sp))
        M, nx = build(Bp, flint.fmpq(0))
        n = len(M)
        A = flint.fmpq_mat(nx, nx, [M[i][j] for i in range(nx) for j in range(nx)])
        W = flint.fmpq_mat(n - nx, nx, [M[i][j] for i in range(nx, n) for j in range(nx)])
        Z = flint.fmpq_mat(n - nx, n - nx, [M[i][j] for i in range(nx, n) for j in range(nx, n)])
        S = A - W.transpose() * Z.solve(W)
        ok1 = all(S[i, j] == Bp.Y(G[i], G[j]) for i in range(nx) for j in range(nx))
        pd = psd_point(M)
        K = Field(list(blk['names']))
        ssym = [K(1)] + list(K.gens) if k == 4 else list(K.gens)
        Bs = Builder(Ring(ssym, K(0), K(1)))
        Ms, _ = build(Bs, K(0))
        Mg = [[(x.n, x.d) for x in row] for row in Ms]
        ob = obstruction(Mg)
        nm, nf = mono_psd(Mg)
        verdict = f'{len(ob[0])} obstructed pairs' if ob else 'diagonal not SF'
        extra = ''
        if ob and ob[0]:
            fk = farkas(Mg, False)
            extra = f'; monomial rescalings: {fk[0]}'
            if fk[0] != 'certificate':
                d2 = rescale_lp(Mg, False)
                if not isinstance(d2, str) and d2 is not None:
                    nm2, nf2 = mono_psd(apply_rescale(Mg, d2, False))
                    extra += f'; after LP rescaling non-PSD {nf2} of {nm2}'
                    nf = min(nf, nf2)
        ''.join(map(str, y))
        print(f'{blk["name"]} column {"".join(map(str, y))}: {nx} states, lift size {n}; Schur complement = column '
              f'group (exact): {ok1}; lift PSD at the point: {pd}; {verdict}; monomial-PSD (q = 1): {nf} non-PSD of '
              f'{nm}{extra}  [{time.time() - t0:.1f}s]', flush=True)
        summary.append((ok1 and pd, nf == 0))
print(f'RESULT k = {k}: {len(summary)} column groups with at least 2 states; lift exact and PSD for '
      f'{sum(1 for a, b in summary if a)}; a fixed-vector representation found for {sum(1 for a, b in summary if b)}')
print('PASS' if all(a for a, b in summary) else 'FAIL')
