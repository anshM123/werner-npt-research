"""The block and the certified objects against ylib.py, the builder from the definitions (M, V_p, D_p, R_p and the
partial transpose, rescaled basis; not Lemma B), exact.

(a) at k = 5, 6 (two rational weight vectors at k = 5, one at k = 6, levels interleaved): 2^{k-2} times the ylib block
    equals l2block.full() entry by entry;
(b) for every size pair (u0, u1) with u0 >= u1 at k = 5, 6 and at (4, 3): the gauge-fixed Schur complements S
    certified by cert_fixed.py (loaded from logs/cert_k*_*_*.pkl), evaluated at q0, equal the projections of the Schur complement of the ylib block
    (exact inverse of its vacuum block) on the test vectors.  Usage: python check_ylib3.py"""
import time
from fractions import Fraction as Fr
import flint
from ylib import Builder, point_ring
from j1lib import to_ylib_state
from r2lib import states_L1
from l2block import Block, ginv
from young import part_data
from saved import load_saved

t0 = time.time()
ok = True


def ylib_block(s):
    k = len(s)
    Bd = Builder(point_ring([flint.fmpq(x.numerator, x.denominator) for x in s]))
    st = states_L1(k)
    ys = [to_ylib_state(x, k) for x in st]
    n = len(st)
    Y = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            v = Bd.Y(ys[i], ys[j])
            Y[i][j] = Y[j][i] = Fr(int(v.p), int(v.q)) * 2 ** (k - 2)
    return Y


for s in ([Fr(1), Fr(4), Fr(1), Fr(4), Fr(4)], [Fr(1, 4), Fr(9), Fr(25, 9), Fr(1), Fr(4)],
          [Fr(1, 9), Fr(1), Fr(1, 9), Fr(1), Fr(1), Fr(1, 9)]):
    Y = ylib_block(s)
    F = Block(s, Fr(1)).full()
    eq = all(Y[i][j] == F[i][j] for i in range(len(F)) for j in range(len(F)))
    ok &= eq
    print(f'(a) k={len(s)} s={[str(x) for x in s]}: 2^(k-2) ylib = l2block: {eq} [{time.time() - t0:.1f}s]', flush=True)


def proj(S, A, B):
    SB = [[sum(S[i][j] * b[j] for j in range(len(b)) if b[j] != 0) for b in B] for i in range(len(S))]
    return [[sum(a[i] * SB[i][l] for i in range(len(a)) if a[i] != 0) for l in range(len(B))] for a in A]


for (u0, u1, q0) in [(5, 0, Fr(1)), (4, 1, Fr(4)), (3, 2, Fr(4)), (3, 2, Fr(1, 9)), (6, 0, Fr(1)), (5, 1, Fr(1, 4)),
                     (4, 2, Fr(9)), (3, 3, Fr(1, 4)), (4, 3, Fr(4))]:
    k = u0 + u1
    s = [Fr(1)] * u0 + [q0] * u1
    Y = ylib_block(s)
    nP = k * (k * (k - 1) // 2)
    YPP = [r[:nP] for r in Y[:nP]]
    YPQ = [r[nP:] for r in Y[:nP]]
    Qi = ginv([r[nP:] for r in Y[nP:]], Fr(1))
    T = [[sum(Qi[a][c] * YPQ[j][c] for c in range(len(Qi))) for j in range(nP)] for a in range(len(Qi))]
    S = [[YPP[i][j] - sum(YPQ[i][a] * T[a][j] for a in range(len(Qi))) for j in range(nP)] for i in range(nP)]
    saved = load_saved(k, u0, u1)
    good = True
    for lam, Ssym in saved.items():
        pd = part_data(k, u0, *lam)
        if pd['W']:
            good &= proj(S, pd['W'], pd['W']) == [[f(q0) for f in row] for row in Ssym]
        if pd['K']:
            good &= all(x == 0 for row in proj(S, pd['K'] + pd['W'], pd['K']) for x in row)
    ok &= good
    print(f'(b) k={k} (u0, u1) = ({u0}, {u1}), q = {q0}: symbolic S of every part = projection of the ylib Schur '
          f'complement, e.K annihilated: {good} [{time.time() - t0:.1f}s]', flush=True)
print('PASS' if ok else 'FAIL')
