"""Lean certificate finder (same mathematics as cert_finder.py, no CVXPY): Clarabel called directly in conic form.

Variables x = (D_1, ..., D_nf, lam).  For each block b (spin L, swap sector, complement of Sym^{2n}) with current
M_b = V diag(e) V^T, w = 1/sqrt(max(|e|, 1)), R = V diag(w):
    S_b(x) = diag(e w^2) + sum_i D_i R^T F_{i,b} R - lam diag(w^2)  in the PSD cone   (exactly M_b(y + D) >= lam I),
    rho -+ D_i >= 0.   Maximise lam.
usage: python cert_lean.py SOLVER LCUT K0 iters n1 n2 ...  (only blocks L <= LCUT enter the SDP; all blocks are checked)"""
import sys, time, pickle, gc
import numpy as np
import scipy.sparse as sp
import clarabel
import scs
from restricted_sdp import setup

def svec_idx(d):
    rows, cols = np.triu_indices(d)          # row-major upper triangle; reorder to column-major
    order = np.lexsort((rows, cols))
    rows, cols = rows[order], cols[order]
    scale = np.where(rows == cols, 1.0, np.sqrt(2.0))
    return rows, cols, scale

def lam_of(data, y):
    lm = np.inf; at = None
    for (L, Q0, Fs) in data:
        M = Q0 + sum(y[i] * Fs[i] for i in range(len(Fs)))
        e = np.linalg.eigvalsh((M + M.T) / 2)[0]
        if e < lm:
            lm, at = e, L
    return lm, at

LCUT = 10**9
def step(data, y, rho):
    data = [blk for blk in data if blk[0] <= LCUT]
    nf = len(y); nx = nf + 1
    Arows = []; bs = []; cones = []
    # box constraints: s = rho - D_i >= 0  ->  A row = +e_i, b = rho ; s = rho + D_i -> A row = -e_i, b = rho
    Abox = sp.vstack([sp.hstack([sp.eye(nf), sp.csc_matrix((nf, 1))]), sp.hstack([-sp.eye(nf), sp.csc_matrix((nf, 1))])])
    Arows.append(Abox); bs.append(np.full(2 * nf, rho)); cones.append(clarabel.NonnegativeConeT(2 * nf))
    for (L, Q0, Fs) in data:
        M = Q0 + sum(y[i] * Fs[i] for i in range(nf))
        M = (M + M.T) / 2
        e, V = np.linalg.eigh(M)
        w = 1.0 / np.sqrt(np.maximum(np.abs(e), 1.0))
        R = V * w[None, :]
        d = len(e)
        r, c, sc = svec_idx(d)
        base = np.zeros((d, d)); np.fill_diagonal(base, e * w * w)
        W = np.zeros((d, d)); np.fill_diagonal(W, w * w)
        cols = []
        for F in Fs:
            Fp = R.T @ F @ R
            Fp = (Fp + Fp.T) / 2
            cols.append(-(Fp[r, c] * sc))           # s = b - A x : A column = -svec(Fp)
        cols.append(W[r, c] * sc)                   # lam enters with -W  ->  A column = +svec(W)
        Arows.append(sp.csc_matrix(np.stack(cols, axis=1)))
        bs.append(base[r, c] * sc)
        cones.append(clarabel.PSDTriangleConeT(d))
    A = sp.vstack(Arows).tocsc()
    b = np.concatenate(bs)
    q = np.zeros(nx); q[-1] = -1.0
    if SOLVER == 'scs':
        # SCS: A x + s = b; our blocks already use s = b - A x with svec of the upper triangle in column-major order,
        # which is the same as the lower triangle in row-major order; SCS wants the lower triangle column-major,
        # i.e. the upper triangle row-major: reorder rows inside each PSD block.
        perm = [np.arange(2 * nf)]
        off = 2 * nf
        sizes = []
        for (L, Q0, Fs) in data:
            d = Q0.shape[0]; sizes.append(d)
            r, c, sc = svec_idx(d)                       # column-major upper
            key = {(i, j): k for k, (i, j) in enumerate(zip(r, c))}
            rr, cc = np.triu_indices(d)                   # row-major upper == column-major lower (transposed)
            perm.append(off + np.array([key[(i, j)] for i, j in zip(rr, cc)]))
            off += len(r)
        perm = np.concatenate(perm)
        A2 = A[perm, :]; b2 = b[perm]
        data_s = dict(A=A2.tocsc(), b=b2, c=q)
        cone = dict(l=2 * nf, s=sizes)
        solver = scs.SCS(data_s, cone, eps_abs=1e-7, eps_rel=1e-7, max_iters=20000, verbose=False)
        sol = solver.solve()
        x = np.array(sol['x']); status = sol['info']['status']
        return x[:nf], x[-1], status
    P = sp.csc_matrix((nx, nx))
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 400
    sol = clarabel.DefaultSolver(P, q, A, b, cones, st).solve()
    x = np.array(sol.x)
    return x[:nf], x[-1], str(sol.status)

import os
RHO0 = float(os.environ.get("RHO0", "4.0"))
def find(n, K0, iters=4, rho=None):
    rho = RHO0 if rho is None else rho
    t0 = time.time()
    free, data = setup(n, K0)
    t1 = time.time()
    y = np.ones(len(free))
    hist = [lam_of(data, y)]
    DMAX = int(os.environ.get("DMAX", "30"))
    def blockmins(yv):
        out = []
        for (L, Q0, Fs) in data:
            M = Q0 + sum(yv[i] * Fs[i] for i in range(len(Fs)))
            out.append(np.linalg.eigvalsh((M + M.T) / 2)[0])
        return np.array(out)
    bm = blockmins(y)
    active = set(i for i, (L, Q0, Fs) in enumerate(data) if Q0.shape[0] <= DMAX or bm[i] < 20)
    for it in range(iters):
        for inner in range(8):
            sub = [data[i] for i in sorted(active)]
            D, lam_pred, status = step(sub, y, rho)
            ynew = y + D
            bm_new = blockmins(ynew)
            viol = [i for i in range(len(data)) if i not in active and bm_new[i] < lam_pred - 1e-6]
            if not viol:
                break
            active |= set(viol)
        lnew = (float(bm_new.min()), data[int(np.argmin(bm_new))][0])
        if lnew[0] > hist[-1][0]:
            y = ynew; hist.append(lnew)
        else:
            rho /= 2; hist.append(hist[-1])
        if it >= 1 and abs(hist[-1][0] - hist[-2][0]) < 1e-4:
            break
    return dict(n=n, K0=K0, free=free, y=y, lam=hist[-1][0], Lat=hist[-1][1], hist=hist, status=status,
                tsetup=t1 - t0, ttotal=time.time() - t0)

SOLVER = 'clarabel'
if __name__ == '__main__':
    SOLVER = sys.argv.pop(1)
    LCUT = int(sys.argv.pop(1))
    K0 = int(sys.argv[1]); iters = int(sys.argv[2])
    for n in [int(x) for x in sys.argv[3:]]:
        r = find(n, K0, iters)
        print(f"n={n} K0={K0}: least eigenvalue off Sym^(2n) = {r['lam']:.5f} (L={r['Lat']}); greedy start "
              f"{r['hist'][0][0]:.4e}; history {[round(float(h[0]), 4) for h in r['hist']]}; {len(r['free'])} free pairs; "
              f"last status {r['status']}; setup {r['tsetup']:.1f}s total {r['ttotal']:.1f}s", flush=True)
        pickle.dump(r, open(f'cert_K{K0}_n{n}.pkl', 'wb'))
        del r; gc.collect()
