"""Symbolic G_z for m = 1 in the gauge x_t = 0 (python-flint rational functions), with a choice of elimination order.

Usage: python closed_t.py u [order] [dehom]
  order: 'qdp' (q, then diagonal X_hh, then off-diagonal X_rh; default), 'qpd', 'old' (row by row, then q).
  dehom: put w_t = 1 (s_t = 1); loses nothing, since G_z(lambda s) = D G_z(s) D with D = diag(1, lambda, ...).
The eliminated block is the (C_u) block on (X_U, q), positive definite on the open orthant (Theorem C), so every pivot
is a nonzero rational function and the result is the Schur complement at every positive weight vector.
Output: closed_t_u{u}[_dehom].pkl with G (list of (num dict, den dict)) and the variable names."""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from closed_gen import build_scaled


def run(u, order='qdp', dehom=False, verbose=True):
    T, U = [0], list(range(1, u + 1))
    t0 = time.time()
    KS, Ms, ti = build_scaled(T, U, dehom)
    if verbose:
        print(f'u = {u}: built and scaled ({time.time() - t0:.1f} s), block size {ti.n}', flush=True)
    tb = [ti.idx[('R', 0)]] + [ti.idx[('Z', 0, h)] for h in U]
    qi = ti.idx[('Q',)]
    diag = [ti.idx[('P', h, h)] for h in U]
    off = [ti.idx[('P', r, h)] for r in U for h in U if r != h]
    if order == 'qdp':
        el = [qi] + diag + off
    elif order == 'qpd':
        el = [qi] + off + diag
    elif order == 'old':
        el = [ti.idx[('P', r, h)] for r in U for h in U] + [qi]
    else:
        raise ValueError(order)
    ordr = el + tb
    A = [[Ms[i][j] for j in ordr] for i in ordr]
    ne, N = len(el), len(ordr)
    t1 = time.time()
    for p in range(ne):
        pv = A[p][p]
        assert pv != 0
        rowp = A[p]
        for i in range(p + 1, N):
            if A[i][p] == 0:
                continue
            f = A[i][p] / pv
            Arow = A[i]
            for j in range(p + 1, N):
                if rowp[j] != 0:
                    Arow[j] = Arow[j] - f * rowp[j]
        if verbose:
            print(f'  pivot {p + 1}/{ne} done ({time.time() - t1:.1f} s)', flush=True)
    nt = len(tb)
    G = [[A[ne + i][ne + j] for j in range(nt)] for i in range(nt)]
    return KS, G, time.time() - t0


if __name__ == '__main__':
    u = int(sys.argv[1])
    order = sys.argv[2] if len(sys.argv) > 2 else 'qdp'
    dehom = len(sys.argv) > 3 and sys.argv[3] == 'dehom'
    KS, G, tt = run(u, order, dehom)
    d = u + 1
    ok = all(G[i][j] == G[j][i] for i in range(d) for j in range(d))
    print('G symmetric:', ok, flush=True)
    names = [f's{i}' for i in range(1 if dehom else 0, u + 1)]
    out = [[(G[i][j].n.to_dict(), G[i][j].d.to_dict()) for j in range(d)] for i in range(d)]
    fn = f'closed_t_u{u}{"_dehom" if dehom else ""}.pkl'
    pickle.dump(dict(G=out, names=names, u=u, m=1, dehom=dehom, order=order), open(fn, 'wb'))
    print(f'total {tt:.1f} s; wrote {fn}')
    print('PASS' if ok else 'FAIL')
