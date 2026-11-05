"""Closed form of G for general (u, m) (arguments), exact (python-flint).  Same method as closed_m1.py.
Optional third argument 'dehom': put w0 = 1 (s0 = 1), which loses nothing since G_z(lambda s) = D G_z(s) D with
D = diag(1, lambda, ..., 1, lambda, ...) positive (the entries are homogeneous of degrees 0, 1, 2).

The type I block of (R_inf)_{k-1} with T = {0}, U = {1, ..., u} (k = u + 1; cf.TypeI, units of 2Y) is built over
Q(w0, ..., wu).  In the scaled coordinates X_{r,h} = w_r w_h x_{r,h}, zeta_{0,h} = w_0 w_h z_h (rho, q unscaled)
every entry is even in every w_i (checked: an odd exponent raises), so the block becomes a matrix over
Q(s0, ..., su), s_i = w_i^2.  The gauge rows x_{1,h} (h in U) are deleted (the kernel vectors p (x) f_h have
P_{1,h}-coordinate w_1 != 0); all other x and q are eliminated by exact Gaussian elimination.  The result is the
Schur complement G (size u + 1) on (rho, z_1, ..., z_u); for m = 1 the whole statement G > 0 is G = G_00 > 0.
Output: closed_m1_u{u}.pkl = dict(G=[[(num dict, den dict)]], names=...)."""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flintfield import Field, F as FF
from cf import TypeI


def halve_poly(p, ctxS):
    d = p.to_dict()
    out = {}
    for m, c in d.items():
        assert all(e % 2 == 0 for e in m), 'odd exponent'
        out[tuple(e // 2 for e in m)] = c
    return ctxS.from_dict(out)


def build_scaled(T, U, dehom=False):
    k = len(T) + len(U)
    if dehom:
        KW = Field([f'w{i}' for i in range(1, k)])
        KS = Field([f's{i}' for i in range(1, k)])
        w = [KW(1)] + list(KW.gens)
    else:
        KW = Field([f'w{i}' for i in range(k)])
        KS = Field([f's{i}' for i in range(k)])
        w = list(KW.gens)
    ti = TypeI(w, T, F=KW)
    M = ti.build()
    n = ti.n
    sc = []
    for lab in ti.labels:
        if lab[0] in ('P', 'Z'):
            sc.append(w[lab[1]] * w[lab[2]])
        else:
            sc.append(KW(1))
    Ms = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            x = M[i][j] * sc[i] * sc[j]
            y = FF(KS, halve_poly(x.n, KS.ctx), halve_poly(x.d, KS.ctx))
            Ms[i][j] = Ms[j][i] = y
    return KS, Ms, ti


def schur(Ms, ti, verbose=True):
    U = ti.U
    r0 = U[0]
    tb = []
    for t in ti.T:
        tb.append(ti.idx[('R', t)])
        for h in U:
            tb.append(ti.idx[('Z', t, h)])
    el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]

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
        if verbose:
            print(f'  pivot {p + 1}/{ne} done ({time.time() - t0:.1f} s)', flush=True)
    nt = len(tb)
    return [[A[ne + i][ne + j] for j in range(nt)] for i in range(nt)]


if __name__ == '__main__':
    u = int(sys.argv[1])
    m = int(sys.argv[2])
    dehom = len(sys.argv) > 3 and sys.argv[3] == 'dehom'
    t0 = time.time()
    T, U = list(range(m)), list(range(m, m + u))
    KS, Ms, ti = build_scaled(T, U, dehom)
    print(f'u = {u}, m = {m}: built and scaled ({time.time() - t0:.1f} s), block size {ti.n}', flush=True)
    G = schur(Ms, ti)
    d = m * (u + 1)
    ok = all(G[i][j] == G[j][i] for i in range(d) for j in range(d))
    print('G symmetric:', ok, flush=True)
    out = [[(G[i][j].n.to_dict(), G[i][j].d.to_dict()) for j in range(d)] for i in range(d)]
    names = [f's{i}' for i in range(1 if dehom else 0, u + m)]
    pickle.dump(dict(G=out, names=names, u=u, m=m, dehom=dehom), open(f'closed_u{u}_m{m}{"_dehom" if dehom else ""}.pkl', 'wb'))
    print(f'total {time.time() - t0:.1f} s')
    print('PASS' if ok else 'FAIL')
