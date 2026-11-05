"""Closed forms of the T-blocks of G for u = 3 (exact, python-flint).

The type I block of (R_inf)_{k-1} with T = {0, 1, 2}, U = {3, 4, 5} (cf.TypeI, units of 2Y) is built over
Q(w0, ..., w5).  In the scaled coordinates X_{r,h} = w_r w_h x_{r,h}, zeta_{t,h} = w_t w_h z_{t,h} (rho, q unscaled)
every entry is even in every w_i (checked), so it is converted to a matrix over Q(s0, ..., s5), s_i = w_i^2.  The
gauge rows x_{3,h} are deleted (the kernel vectors p (x) f_h have P_{3,h}-coordinate w_3 != 0), and all other x and q
are eliminated by exact Gaussian elimination over Q(s0, ..., s5).  The result is the Schur complement G on
(rho_t, z_{t,3}, z_{t,4}, z_{t,5}), t = 0, 1, 2.  Output (closed_u3.pkl): the 4 x 4 blocks G00, G01 as (numerator,
denominator) dicts over (s0, ..., s5); also G11, G22, G02, G12 for the symmetry checks."""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FF
from cf import TypeI

NAMES_W = ['w0', 'w1', 'w2', 'w3', 'w4', 'w5']
NAMES_S = ['s0', 's1', 's2', 's3', 's4', 's5']


def halve_poly(p, ctxS):
    d = p.to_dict()
    out = {}
    for m, c in d.items():
        assert all(e % 2 == 0 for e in m), 'odd exponent'
        out[tuple(e // 2 for e in m)] = c
    return ctxS.from_dict(out) if out else ctxS.from_dict({})


def build_scaled(T, U):
    """Return (KS, Ms, labels, keep) with Ms the scaled block over Q(s)."""
    k = len(T) + len(U)
    KW = Field(NAMES_W[:k])
    KS = Field(NAMES_S[:k])
    w = list(KW.gens)
    ti = TypeI(w, T, F=KW)
    M = ti.build()
    n = ti.n
    # scaling vector
    sc = []
    for lab in ti.labels:
        if lab[0] == 'P':
            sc.append(w[lab[1]] * w[lab[2]])
        elif lab[0] == 'Z':
            sc.append(w[lab[1]] * w[lab[2]])
        else:
            sc.append(KW(1))
    Ms = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            x = M[i][j] * sc[i] * sc[j]
            num = halve_poly(x.n, KS.ctx)
            den = halve_poly(x.d, KS.ctx)
            y = FF(KS, num, den)
            Ms[i][j] = Ms[j][i] = y
    return KS, Ms, ti


def schur(KS, Ms, ti):
    U = ti.U
    r0 = U[0]
    tb = []
    for t in ti.T:
        tb.append(ti.idx[('R', t)])
        for h in U:
            tb.append(ti.idx[('Z', t, h)])
    el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
    # eliminate X_T first, then X_U, then q
    def key(i):
        lab = ti.labels[i]
        if lab[0] == 'P' and lab[1] in ti.T:
            return (0, i)
        if lab[0] == 'P':
            return (1, i)
        return (2, i)
    el.sort(key=key)
    order = el + tb
    A = [[Ms[i][j] for j in order] for i in order]
    ne = len(el)
    N = len(order)
    t0 = time.time()
    for p in range(ne):
        pv = A[p][p]
        assert pv != 0
        rowp = A[p]
        for i in range(p + 1, N):
            if A[i][p] == 0:
                continue
            f = A[i][p] / pv
            Ai = A[i]
            for j in range(p + 1, N):
                if rowp[j] != 0:
                    Ai[j] = Ai[j] - f * rowp[j]
        print(f'  pivot {p + 1}/{ne} done ({time.time() - t0:.1f} s)', flush=True)
    nt = len(tb)
    return [[A[ne + i][ne + j] for j in range(nt)] for i in range(nt)]


if __name__ == '__main__':
    t0 = time.time()
    T, U = [0, 1, 2], [3, 4, 5]
    KS, Ms, ti = build_scaled(T, U)
    print(f'built and scaled ({time.time() - t0:.1f} s), size {ti.n}', flush=True)
    G = schur(KS, Ms, ti)
    print(f'Schur complement done ({time.time() - t0:.1f} s)', flush=True)
    d = 4
    ok = all(G[i][j] == G[j][i] for i in range(12) for j in range(12))
    print('G symmetric:', ok, flush=True)
    out = {}
    for name, (bi, bj) in (('G00', (0, 0)), ('G01', (0, 1)), ('G11', (1, 1)), ('G22', (2, 2)), ('G02', (0, 2)), ('G12', (1, 2))):
        out[name] = [[(G[d * bi + i][d * bj + j].n.to_dict(), G[d * bi + i][d * bj + j].d.to_dict()) for j in range(d)] for i in range(d)]
    pickle.dump(out, open('closed_u3.pkl', 'wb'))
    for name in ('G00', 'G01'):
        for i in range(d):
            for j in range(d):
                n_, d_ = out[name][i][j]
                print(name, i, j, 'num terms', len(n_), 'den terms', len(d_),
                      'num deg', max(sum(m) for m in n_) if n_ else None, 'den deg', max(sum(m) for m in d_), flush=True)
    print(f'total {time.time() - t0:.1f} s')
    print('PASS' if ok else 'FAIL')
