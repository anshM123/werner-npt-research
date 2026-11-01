"""Pilot exact certificate: q_2(P) >= 0 for every rank-2 projector P on (C^3)^(x)2 (t = +1), via the two-letter
four-replica relaxation (X >= 0, moments, swap identities built in, ORTH, MARG, covariance blocks).

Pipeline
 1. exact data (exact_core.ExactRelax2): variables x in Q^55, objective c, moment rows (E_m x = e_m), identity rows
    (E_i x = 0), PSD blocks M_k(x) (rational symmetric matrices, linear in x);
 2. exact zeros: u = a (x) psi, v = a' (x) psi with rational orthonormal a, a' (Cayley transforms of integer skew
    matrices) and integer psi, site-symmetrised; x* = sum of K such points; exact kernels N_k of M_k(x*);
 3. floating-point solve of the face-reduced problem  min c.x  s.t.  E x = e,  N_k^T M_k(x) N_k >= 0  (Clarabel);
 4. rounding: the dual (y, S_k) is rounded to rationals and corrected exactly (least-norm correction) so that
        c = E^T y + sum_k M_k^*(N_k S_k N_k^T)   and   e . y = 0
    hold exactly (all 55 + 1 equations);
 5. exact LDL^T check that every S_k is positive semidefinite.
If all checks pass, for every rank-2 projector P (true point x(P) of the relaxation):
    q_2(P) = c . x(P) = e . y + sum_k Tr(S_k N_k^T M_k(x(P)) N_k) >= 0.
usage: python exact_pilot_n2.py [K]"""
import sys, time, itertools, numpy as np, flint, scipy.sparse as sp, clarabel
from fractions import Fraction
from exact_core import *

def log(s):
    print(s, flush=True)

def cayley(rng):
    """rational orthogonal 3x3 matrix (I - S)(I + S)^{-1}, S integer skew-symmetric"""
    a, b, c = [int(z) for z in rng.integers(-3, 4, size=3)]
    S = flint.fmpq_mat(3, 3, [0, a, b, -a, 0, c, -b, -c, 0])
    I = eye(3)
    return (I - S) * (I + S).inv()

def exact_psd(S):
    """exact LDL^T test of a rational symmetric matrix (list of lists of fmpq): returns (is_psd, smallest pivot)"""
    k = len(S); M = [row[:] for row in S]; minpiv = None
    for i in range(k):
        p = M[i][i]
        if p < 0: return False, p
        if p == 0:
            if any(M[i][j] != 0 for j in range(i + 1, k)): return False, p
            continue
        minpiv = p if minpiv is None or p < minpiv else minpiv
        for j in range(i + 1, k):
            if M[j][i] != 0:
                f = M[j][i] / p
                for l in range(i + 1, k): M[j][l] -= f * M[i][l]
    return True, minpiv

def main(K=4):
    t0 = time.time()
    A, B1, A3 = site_algebras()
    rel = ExactRelax2(A, B1, A3, log=lambda s: None)
    N = rel.N
    log(f"[1] exact data: {N} parameters, {rel.ident.nrows()} identity rows, {len(rel.psd)} PSD blocks; {time.time() - t0:.1f}s")
    # ---- moments: independent rows (exact)
    Em, em = rel.mom
    aug = vstack([hstack([Em, em]), hstack([rel.ident, zeros(rel.ident.nrows(), 1)])])
    R, rk = aug.rref()
    E = flint.fmpq_mat(rk, N, [R[i, j] for i in range(rk) for j in range(N)])
    e = flint.fmpq_mat(rk, 1, [R[i, N] for i in range(rk)])
    # consistency: no row of the form 0 = nonzero
    assert all(any(E[i, j] != 0 for j in range(N)) for i in range(rk)), "inconsistent equations"
    log(f"    equalities: {Em.nrows()} moment rows and {rel.ident.nrows()} identity rows, {rk} independent")
    # ---- exact zeros and kernels
    rng = np.random.default_rng(2026)
    xs = []
    for k in range(K):
        Qm = cayley(rng); a = [Qm[0, j] for j in range(3)]; a2 = [Qm[1, j] for j in range(3)]
        psi = [fq(int(z)) for z in rng.integers(-4, 5, size=3)]
        if all(z == 0 for z in psi): psi[0] = fq(1)
        hz = {}
        for c in FAMS:
            h1 = rel.site_hat([a if ch == 'u' else a2 for ch in c]); h2 = rel.site_hat([psi] * 4)
            hz[c] = {(l1, l2): (kron(h1[l1], h2[l2]) + kron(h2[l1], h1[l2])) * fq(1, 2) for l1 in range(4) for l2 in range(4)}
        x, ok = rel.x_from_hats(hz)
        mres, ires, obj, mins = rel.evaluate(x)
        assert ok and obj == 0 and all(mres[i, 0] == 0 for i in range(mres.nrows())) and all(ires[i, 0] == 0 for i in range(ires.nrows())), "bad zero"
        xs.append(x)
    xstar = xs[0]
    for x in xs[1:]: xstar = xstar + x
    blocks = []      # (name, size, Nk (exact), list over params of Nk^T M_k(e_i) Nk (exact))
    dims = [0, 0]
    for name, size, mats in rel.psd:
        Bx = zeros(size, size)
        for i in range(N):
            if xstar[i, 0] != 0: Bx = Bx + mats[i] * xstar[i, 0]
        Nk = nullspace_q(Bx)
        dims[0] += size; dims[1] += Nk.ncols()
        if Nk.ncols() == 0: continue
        red = [Nk.transpose() * mats[i] * Nk for i in range(N)]
        blocks.append((name, Nk.ncols(), Nk, red))
    log(f"[2] {K} exact zeros (objective 0, moments and identities exact); face: block dimension {dims[0]} -> {dims[1]}, "
        f"{len(blocks)} blocks kept; {time.time() - t0:.1f}s")
    # ---- floating-point solve of the reduced problem
    Ef = to_np(E); ef = to_np(e).ravel(); cf = to_np(rel.c).ravel()
    rows = [Ef]; b = list(ef); cones = [clarabel.ZeroConeT(Ef.shape[0])]
    from rep4 import svec, unsvec, svec_pairs
    nnrows = []; psd_rows = []
    for name, k, Nk, red in blocks:
        Mf = np.array([svec(to_np(r)) for r in red]).T        # svec rows x N
        if k == 1: nnrows.append(-Mf[0])
        else: psd_rows.append((k, -Mf))
    if nnrows:
        rows.append(np.array(nnrows)); b += [0.0] * len(nnrows); cones.append(clarabel.NonnegativeConeT(len(nnrows)))
    for k, Mf in psd_rows:
        rows.append(Mf); b += [0.0] * Mf.shape[0]; cones.append(clarabel.PSDTriangleConeT(k))
    Amat = sp.csc_matrix(np.vstack(rows)); bvec = np.array(b)
    st = clarabel.DefaultSettings(); st.verbose = False; st.max_iter = 300
    st.tol_gap_abs = 1e-11; st.tol_gap_rel = 1e-11; st.tol_feas = 1e-11
    sol = clarabel.DefaultSolver(sp.csc_matrix((N, N)), cf, Amat, bvec, cones, st).solve()
    z = np.array(sol.z)
    log(f"[3] float solve: status {sol.status}, primal {sol.obj_val:.3e}, dual {-(bvec @ z):.3e}, "
        f"||A^T z + c|| {np.abs(Amat.T @ z + cf).max():.1e}; {time.time() - t0:.1f}s")
    # dual variables: y = -z_eq, S_k from the cone slices
    neq = Ef.shape[0]; y0 = -z[:neq]; row = neq; S0 = []
    if nnrows:
        for j in range(len(nnrows)): S0.append(np.array([[z[row + j]]]))
        row += len(nnrows)
    for k, Mf in psd_rows:
        kk = Mf.shape[0]; S0.append(unsvec(z[row:row + kk], k)); row += kk
    order = [blk for blk in blocks if blk[1] == 1] + [blk for blk in blocks if blk[1] > 1]
    mineig = min(np.linalg.eigvalsh(S).min() for S in S0)
    log(f"    smallest eigenvalue of the reduced dual blocks: {mineig:.2e}")
    if mineig <= 1e-9:
        log("FAILED: the reduced dual is not strictly feasible (face not minimal); rounding not attempted"); return
    # ---- rounding and exact least-norm correction
    SCALE = 2 ** 40
    rq = lambda v: fq(int(round(float(v) * SCALE)), SCALE)
    unk = []          # list of (kind, index) with kind 'y' or ('S', block, a, b)
    for i in range(neq): unk.append(('y', i))
    for bi, (name, k, Nk, red) in enumerate(order):
        for a_ in range(k):
            for b_ in range(a_, k): unk.append(('S', bi, a_, b_))
    nu = len(unk)
    # equation matrix: rows 0..N-1: sum_e y_e E[e,i] + sum_k Tr(S_k red_k[i]) = c_i ;  row N: sum_e y_e e_e = 0
    Aeq = flint.fmpq_mat(N + 1, nu); rhs = flint.fmpq_mat(N + 1, 1)
    for i in range(N): rhs[i, 0] = rel.c[0, i]
    z0 = flint.fmpq_mat(nu, 1)
    for u_, kd in enumerate(unk):
        if kd[0] == 'y':
            ei = kd[1]
            for i in range(N): Aeq[i, u_] = E[ei, i]
            Aeq[N, u_] = e[ei, 0]
            z0[u_, 0] = rq(y0[ei])
        else:
            _, bi, a_, b_ = kd; red = order[bi][3]; w = 1 if a_ == b_ else 2
            for i in range(N): Aeq[i, u_] = red[i][a_, b_] * w
            z0[u_, 0] = rq(S0[bi][a_, b_])
    resid = rhs - Aeq * z0
    rmax = max(abs(float(resid[i, 0])) for i in range(N + 1))
    # least-norm correction: z = z0 + Aeq^T (Aeq Aeq^T)^{-1} resid, on an independent set of rows
    R, rk = Aeq.transpose().rref()       # independent equations = pivot columns of Aeq^T
    piv = []; r_ = 0
    for j in range(Aeq.nrows()):
        if r_ < rk and R[r_, j] == 1 and all(R[ii, j] == 0 for ii in range(rk) if ii != r_): piv.append(j); r_ += 1
    Ai = flint.fmpq_mat(len(piv), nu, [Aeq[i, j] for i in piv for j in range(nu)])
    ri = flint.fmpq_mat(len(piv), 1, [resid[i, 0] for i in piv])
    corr = Ai.transpose() * ((Ai * Ai.transpose()).solve(ri))
    zc = z0 + corr
    final = rhs - Aeq * zc
    exact_ok = all(final[i, 0] == 0 for i in range(N + 1))
    cmax = max(abs(float(corr[i, 0])) for i in range(nu))
    log(f"[4] rounding (2^-40 grid): residual before correction {rmax:.1e}; rank of the equations {rk} of {N + 1}; "
        f"correction size {cmax:.1e}; identity exact after correction: {exact_ok}; {time.time() - t0:.1f}s")
    if not exact_ok:
        log("FAILED: the corrected dual does not satisfy the identity exactly (dependent equations inconsistent)"); return
    # ---- exact PSD checks
    Sq = []; off = neq
    for bi, (name, k, Nk, red) in enumerate(order):
        S = [[fq(0)] * k for _ in range(k)]
        for a_ in range(k):
            for b_ in range(a_, k):
                S[a_][b_] = zc[off, 0]; S[b_][a_] = zc[off, 0]; off += 1
        Sq.append((name, S))
    allpsd = True; minpiv = None
    for name, S in Sq:
        ok, mp = exact_psd(S)
        if not ok:
            allpsd = False; log(f"    block {name}: NOT PSD (pivot {mp})"); break
        if mp is not None and (minpiv is None or mp < minpiv): minpiv = mp
    dual_obj = sum((e[i, 0] * zc[i, 0] for i in range(neq)), fq(0))
    # sizes of the rationals in the certificate
    nb = 0; db = 0
    for i in range(nu):
        qv = zc[i, 0]; nb = max(nb, int(abs(qv.p)).bit_length()); db = max(db, int(qv.q).bit_length())
    nbd, dbd = bits(Aeq)
    log(f"[5] exact LDL^T: all {len(Sq)} reduced dual blocks PSD: {allpsd} (smallest pivot {float(minpiv) if minpiv is not None else 0:.2e}); "
        f"dual objective e.y = {dual_obj}; {time.time() - t0:.1f}s")
    log(f"    rational sizes: certificate entries up to {nb} numerator bits / {db} denominator bits; "
        f"data entries up to {nbd} / {dbd} bits")
    # independent check on an exact entangled true point: q_2(P) = e.y + sum_k Tr(S_k N_k^T M_k(x(P)) N_k), terms >= 0
    rng2 = np.random.default_rng(7)
    p = rng2.integers(-3, 4, size=9); q = rng2.integers(-3, 4, size=9)
    u = p; v = (p @ p) * q - (p @ q) * p
    L = {'u': u.reshape(3, 3), 'v': v.reshape(3, 3)}
    hats = {c: two_site_hats(rel, [L[ch] for ch in c]) for c in FAMS}
    xP, okP = rel.x_from_hats(hats)
    Pm = np.outer(u, u).astype(object) * Fraction(1, int(u @ u)) + np.outer(v, v).astype(object) * Fraction(1, int(v @ v))
    q2 = exact_q2(Pm)
    lhs = (rel.c * xP)[0, 0]
    terms = []
    for (name, S), (nm2, k, Nk, red) in zip(Sq, order):
        Mx = zeros(k, k)
        for i in range(N):
            if xP[i, 0] != 0: Mx = Mx + red[i] * xP[i, 0]
        terms.append(sum((S[a_][b_] * Mx[b_, a_] for a_ in range(k) for b_ in range(k)), fq(0)))
    eyx = sum(((E * xP)[i, 0] * zc[i, 0] for i in range(neq)), fq(0))
    rhs_val = eyx + sum(terms, fq(0))
    allpos = all(tm >= 0 for tm in terms)
    log(f"[6] check on an entangled true point: c.x(P) = {lhs} = q_2(P) = {q2}: {Fraction(int(lhs.p), int(lhs.q)) == q2 and okP}; "
        f"y.E x(P) = {eyx}; c.x(P) == y.E x(P) + sum_k Tr(S_k M_k): {lhs == rhs_val}; all {len(terms)} terms >= 0: {allpos} "
        f"(smallest {float(min(terms)):.2e})")
    if allpsd and exact_ok and dual_obj == 0 and lhs == rhs_val and allpos:
        print(f"n=2 t=+1 pilot: c = E^T y + sum_k M_k^*(N_k S_k N_k^T) exactly, e.y = 0, all S_k PSD (exact LDL^T), "
              f"{len(Sq)} blocks, max {nb}/{db} bits -> CERTIFIED", flush=True)
    else:
        print(f"n=2 t=+1 pilot: FAILED (psd {allpsd}, identity {exact_ok}, dual objective {dual_obj})", flush=True)

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
