"""Shared helpers: loading the pickled closed forms of equal_m.py into python-flint, exact evaluation, and the exact
Schur complement G (z-coordinates) of a whole type I block of (R_inf)_{k-1} for any T, U at rational weights."""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from cf import TypeI

Q = flint.fmpq
NAMES = ('u', 'sigma', 's', 's1', 's2')
CTX = flint.fmpq_mpoly_ctx.get(NAMES, 'lex')


def poly(d, ctx=CTX):
    return ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})


def load(path='equal_m.pkl'):
    R = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), path), 'rb'))
    out = {}
    for key in ('inv', 'std'):
        G = R[key]['G']
        out[key] = [[(poly(n), poly(d)) for (n, d) in row] for row in G]
    return out


def ev(nd, vals):
    n, d = nd
    vals = [Q(v) for v in vals]
    return n(*vals) / d(*vals)


def schur_G(wT, wU, gauge_r0=None):
    """Exact G in the coordinates (rho_t, zeta_{t,h}) (order: per t, rho_t then zeta_{t,h}, h in U), T = copies
    0..m-1 with amplitudes wT, U = the next u copies with amplitudes wU; gauge: the coordinates P_{r0,h} deleted
    (r0 = first copy of U), a complement of the kernel span{p (x) f_h}.  Returns (G, labels)."""
    w = [Q(x) for x in list(wT) + list(wU)]
    m = len(wT)
    ti = TypeI(w, list(range(m)), F=Q)
    M = ti.build()
    U = ti.U
    r0 = U[0] if gauge_r0 is None else gauge_r0
    tb = []
    for t in range(m):
        tb.append(ti.idx[('R', t)])
        for h in U:
            tb.append(ti.idx[('Z', t, h)])
    el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
    C = flint.fmpq_mat([[M[i][j] for j in el] for i in el])
    B = flint.fmpq_mat([[M[i][j] for j in tb] for i in el])
    A = flint.fmpq_mat([[M[i][j] for j in tb] for i in tb])
    G = A - B.transpose() * C.solve(B)
    return G, w, m, U


def Gz(wT, wU, gauge_r0=None):
    """z-coordinates: zeta_{t,h} = w_t w_h z_{t,h}."""
    G, w, m, U = schur_G(wT, wU, gauge_r0)
    u = len(U)
    d = []
    for t in range(m):
        d.append(Q(1))
        for h in U:
            d.append(w[t] * w[h])
    n = len(d)
    return flint.fmpq_mat([[G[i, j] * d[i] * d[j] for j in range(n)] for i in range(n)]), m, u


def ldl_pivots(Mq):
    n = Mq.nrows()
    A = [[Mq[i, j] for j in range(n)] for i in range(n)]
    piv = []
    for p in range(n):
        pv = A[p][p]
        piv.append(pv)
        if pv == 0:
            return piv
        for i in range(p + 1, n):
            f = A[i][p] / pv
            if f != 0:
                for j in range(p + 1, n):
                    A[i][j] -= f * A[p][j]
    return piv


# ---------------------------------------------------------------------------------------------------------------
# Direct exact computations from the Lemma B builder (fastj1.FastJ1; independent of cf.py and of equal_m.py)
def lemmab_block(w, T, bscale=1):
    """whole type I block of (R_inf)_{k-1} for the set T, units 2Y, from the Lemma B builder; states in the order of
    fastdirect.typeI_states (= the order of cf.TypeI): P_{r,h}, R_t, Z_{t,h}, Q.  bscale multiplies the B-part."""
    from fractions import Fraction as Fr
    from fastj1 import FastJ1
    from fastdirect import typeI_states
    k = len(w)
    st, lab = typeI_states(k, T)
    conv = []
    for row, col in st:
        nz = [c for c in range(k) if row[c] != 0]
        conv.append((None if not nz else (row[nz[0]], nz[0]), col))
    M = FastJ1([Fr(x) for x in w], bscale=bscale).matrix(conv)
    return M, lab


def qf(x):
    from fractions import Fraction as Fr
    return Q(x.numerator, x.denominator) if isinstance(x, Fr) else Q(x)


def direct_G(wT, om, u, bscale=1):
    """T = copies 0..m-1 (amplitudes wT), U = u copies of amplitude om.  Returns (Gz, pivots of the gauge-deleted
    block, the block, labels).  Gz: z-coordinates zeta_{t,h} = w_t om z_{t,h}, order per t: rho_t, z_{t,h}."""
    from fractions import Fraction as Fr
    m = len(wT)
    w = [Fr(x) for x in wT] + [Fr(om)] * u
    k = len(w)
    T = list(range(m))
    U = list(range(m, k))
    Mb, lab = lemmab_block(w, T, bscale)
    n = len(lab)
    idx = {l: i for i, l in enumerate(lab)}
    r0 = U[0]
    keep = [i for i in range(n) if not (lab[i][0] == 'P' and lab[i][1] == r0)]
    piv = ldl_pivots(flint.fmpq_mat([[qf(Mb[i][j]) for j in keep] for i in keep]))
    tb = []
    for t in T:
        tb.append(idx[('R', t)])
        for h in U:
            tb.append(idx[('Z', t, h)])
    el = [i for i in keep if i not in tb]
    C = flint.fmpq_mat([[qf(Mb[i][j]) for j in el] for i in el])
    B = flint.fmpq_mat([[qf(Mb[i][j]) for j in tb] for i in el])
    A = flint.fmpq_mat([[qf(Mb[i][j]) for j in tb] for i in tb])
    G = A - B.transpose() * C.solve(B)
    d = []
    for t in T:
        d.append(Q(1))
        for h in U:
            d.append(qf(Fr(wT[t])) * qf(Fr(om)))
    Gz = [[G[i, j] * d[i] * d[j] for j in range(len(d))] for i in range(len(d))]
    return Gz, piv, Mb, lab


def direct_conditions(wT, om, u, H, Gz=None):
    """the five Gershgorin quantities of prove_em.py computed from the exact G of the whole block (no closed form):
    dict name -> list of values (per copy t for B, per ordered pair (t, t') for A)."""
    from fractions import Fraction as Fr
    if Gz is None:
        Gz = direct_G(wT, om, u)[0]
    m = len(wT)
    D = u + 1
    s = qf(Fr(om)) ** 2
    sT = [qf(Fr(x)) ** 2 for x in wT]
    sig = sum(sT, Q(0))
    S = u * s
    H = Q(H)

    def inv_blk(t, tp):
        b0, b1 = t * D, tp * D
        return [[Gz[b0][b1], sum((Gz[b0][b1 + 1 + h] for h in range(u)), Q(0))],
                [sum((Gz[b0 + 1 + h][b1] for h in range(u)), Q(0)),
                 sum((Gz[b0 + 1 + h][b1 + 1 + hp] for h in range(u) for hp in range(u)), Q(0))]]

    def std_blk(t, tp):
        b0, b1 = t * D, tp * D
        return Gz[b0 + 1][b1 + 1] - Gz[b0 + 1][b1 + 2]
    out = dict(A_inv=[], A_std=[], B_inv_00=[], B_inv_det=[], B_std=[])
    kap2 = [st * S * sig / (sig + S) for st in sT]
    for t in range(m):
        th1 = 1 - sT[t] / sig
        g = inv_blk(t, t)
        out['B_inv_00'].append(g[0][0] - H * th1)
        out['B_inv_det'].append((g[0][0] - H * th1) * (g[1][1] - H * th1 * kap2[t]) - g[0][1] * g[1][0])
        out['B_std'].append(std_blk(t, t) / (sT[t] * s) - H * th1)
        for tp in range(m):
            if tp == t:
                continue
            g = inv_blk(t, tp)
            fro = g[0][0] ** 2 + g[0][1] ** 2 / kap2[tp] + g[1][0] ** 2 / kap2[t] + g[1][1] ** 2 / (kap2[t] * kap2[tp])
            tgt = H * H * sT[t] * sT[tp] / sig ** 2
            out['A_inv'].append(tgt - fro)
            out['A_std'].append(tgt - std_blk(t, tp) ** 2 / (sT[t] * sT[tp] * s * s))
    return out
