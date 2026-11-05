"""The reduced forms of l1sym.py as identities of rational functions in the weight ratio, at fixed level sizes, against
the independent builder ylib.Builder (which builds Y = A - B from the definitions: M, V_p, D, R and the partial
transpose; it does not use Lemma B), over Q(q) with squared weights 1 on level 0 and q on level 1.

For (u0, u1) in {(2,2), (3,1), (1,3), (3,2), (2,3), (4,1), (1,4)} (k = 4, 5): the letter-1 block of (R_inf)_1 is built by
ylib over Q(q) (in units 2^{k-1}: ylib's rescaled block times 2^{k-1}); every part that exists at these sizes (irrep on
level 0) is projected on the explicit test vectors; R_PP, R_PQ, R_QQ and the Schur complement S of the vacuum states
(computed from the whole block, then projected) must equal the forms of l1sym.py computed over Q(q) at the same fixed
sizes, identically in q.  The block's kernel vectors (indicators 1_c) are checked identically in q as well.
Usage: python check_ylib_sym.py"""
import time
from fractions import Fraction as Fr
import flint
from flintfield import Field
from ylib import Builder, Ring
from j1lib import block_states, to_ylib_state
from l1sym import Model, reduce_part, keep_for, part_spec, mm, mT, minv
from check_reduce import concrete_labels

t0 = time.time()
K = Field(['q'])
q = K.gens[0]
one, zero = K(1), K(0)
ok = True
for (u0, u1) in [(2, 2), (3, 1), (1, 3), (3, 2), (2, 3), (4, 1), (1, 4)]:
    k = u0 + u1
    s = [one] * u0 + [q] * u1
    Bd = Builder(Ring(s, zero, one))
    st = block_states(k, ())
    ys = [to_ylib_state(x, k) for x in st]
    n = len(st)
    sc = K(2 ** (k - 1))
    Mt = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            Mt[i][j] = Mt[j][i] = Bd.Y(ys[i], ys[j]) * sc
    sym = all(Bd.Y(ys[j], ys[i]) * sc == Mt[i][j] for i in range(n) for j in range(i + 1, min(n, i + 4)))
    Pi = [i for i, (a, y) in enumerate(st) if a is not None]
    Qi = [i for i, (a, y) in enumerate(st) if a is None]
    kerok = True
    for c in range(k):
        v = [1 if (a is not None and y[c] == 1) else 0 for (a, y) in st]
        kerok &= all(sum((Mt[i][j] for j in range(n) if v[j]), zero) == 0 for i in range(n))
    A = [[Mt[i][j] for j in Qi] for i in Qi]
    Ai = minv(A, one)
    BPQ = [[Mt[i][j] for j in Qi] for i in Pi]
    C = mm(mm(BPQ, Ai, zero), mT(BPQ), zero)
    SP = [[Mt[Pi[a]][Pi[c]] - C[a][c] for c in range(len(Pi))] for a in range(len(Pi))]
    print(f'(u0, u1) = ({u0}, {u1}): ylib block over Q(q) of size {n}, symmetric {sym}, kernel vectors 1_c '
          f'identically {kerok} [{time.time() - t0:.1f}s]', flush=True)
    ok &= sym and kerok
    for part in ('inv', 'std', 'two', 'anti', 'mix'):
        if (part == 'std' and u0 < 2) or (part == 'two' and u0 < 4) or (part == 'anti' and u0 < 3) or \
                (part == 'mix' and (u0 < 2 or u1 < 2)):
            continue
        kp, kq, gauge = keep_for(part, u0, u1)
        red = reduce_part(Model(one, [one, q], [K(u0), K(u1)]), part, kp, kq)
        labs = concrete_labels(part, u0, u1)
        named, P, Qv, ker, names = part_spec(part)
        P = [P[i] for i in kp]
        Qv = [Qv[i] for i in kq] if (Qv and kq is not None) else []
        VP = [[f(labs[st[i][0][1]], labs[st[i][1].index(1)]) if i in Pi else 0 for i in range(n)] for f in P]
        VQ = []
        for g in Qv:
            v = [0] * n
            for i in Qi:
                pr = [t for t in range(k) if st[i][1][t] == 1]
                v[i] = g(labs[pr[0]], labs[pr[1]])
            VQ.append(v)

        def form(X, Y):
            return [[sum((Mt[i][j] * (x[i] * y[j]) for i in range(n) if x[i] for j in range(n) if y[j]), zero)
                     for y in Y] for x in X]
        same = form(VP, VP) == red['RPP']
        if VQ:
            same &= form(VP, VQ) == red['RPQ'] and form(VQ, VQ) == red['RQQ']
        VPr = [[v[i] for i in Pi] for v in VP]
        Sx = [[sum((SP[a][c] * (x[a] * y[c]) for a in range(len(Pi)) if x[a] for c in range(len(Pi)) if y[c]), zero)
               for y in VPr] for x in VPr]
        sameS = Sx == red['S']
        ok &= same and sameS
        print(f'  {part:4s}: R_PP, R_PQ, R_QQ equal identically in q {same}; S equal identically in q {sameS} '
              f'[{time.time() - t0:.1f}s]', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
