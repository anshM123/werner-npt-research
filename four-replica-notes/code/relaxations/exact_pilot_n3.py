"""Exact certificate attempt at n = 3 (t = +1, rank-2 projectors) with the two-letter relaxation
(X >= 0, moments, swap identities built in, ORTH, MARG, covariance blocks), same pipeline as exact_pilot_n2.py with:
  * n-site exact data (exact_coreN.ExactRelaxN), exact word moments of integer vectors via primes + CRT;
  * exact zeros: k = 1 (a (x) psi, a' (x) psi), k = 3 (orthogonal products), and the interpolating family
    a (x) (3 w+ + 4 s w-), a' (x) (4 w+ + 3 s w-) followed by a rational rotation of the (u, v) basis;
  * rounding with an exact correction on a column subset of full row rank (pivoted QR on the float matrix).
usage: python exact_pilot_n3.py"""
import sys, time, itertools, numpy as np, flint, scipy.sparse as sp, scipy.linalg, clarabel
from gmpy2 import mpq
from exact_core import site_algebras, FAMS
from exact_coreN import *
from exact_site import nullspace_q
from rep4 import svec, unsvec

def log(s): print(s, flush=True)

def cayley_int(rng):
    """integer matrix d * Q with Q rational orthogonal (Cayley transform of an integer skew matrix); returns (M, d)"""
    a, b, c = [int(z) for z in rng.integers(-2, 3, size=3)]
    S = flint.fmpq_mat(3, 3, [0, a, b, -a, 0, c, -b, -c, 0]); I = flint.fmpq_mat(3, 3, [1, 0, 0, 0, 1, 0, 0, 0, 1])
    Qm = (I - S) * (I + S).inv()
    d = 1
    for i in range(3):
        for j in range(3): d = d * int(Qm[i, j].q) // np.gcd(d, int(Qm[i, j].q))
    M = np.array([[int(Qm[i, j] * d) for j in range(3)] for i in range(3)], dtype=object)
    return M, d

def zero_vectors(kind, rng):
    if kind == 'k1':
        M, d = cayley_int(rng); psi = rng.integers(-3, 4, size=9).astype(object)
        if not any(psi): psi[0] = 1
        return np.kron(M[0], psi), np.kron(M[1], psi)
    if kind == 'k3':
        u = np.array([1], dtype=object); v = np.array([1], dtype=object)
        for s in range(3):
            M, d = cayley_int(rng); u = np.kron(u, M[0]); v = np.kron(v, M[1])
        return u, v
    # interpolating + rotation
    Ma, _ = cayley_int(rng); Mb, _ = cayley_int(rng); Mc, _ = cayley_int(rng)
    wp = np.kron(Mb[0], Mc[0]); wm = np.kron(Mb[1], Mc[1]); s = 1 if rng.integers(0, 2) else -1
    phi = 3 * wp + 4 * s * wm; chi = 4 * wp + 3 * s * wm
    u = np.kron(Ma[0], phi); v = np.kron(Ma[1], chi)
    r1, r2 = [(3, 4), (5, 12), (8, 15)][int(rng.integers(0, 3))]
    return r1 * u + r2 * v, -r2 * u + r1 * v

def exact_psd(S):
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

def main():
    t0 = time.time()
    A, B1, A3 = site_algebras()
    rel = ExactRelaxN(3, A, B1, A3, log=lambda s: None)
    rel.build(log=lambda s: None)
    N = rel.N
    log(f"[1] exact data n=3: {N} parameters, {rel.ident.shape[0]} identity rows, {len(rel.psd)} PSD blocks; {time.time() - t0:.0f}s")
    # independent equalities (moments + identities), exact
    Mm, em = rel.mom
    aug = np.vstack([np.hstack([Mm, em.reshape(-1, 1)]), np.hstack([rel.ident, zeros_o((rel.ident.shape[0], 1))])])
    Rr, rk = np2fm(aug).rref()
    E = np.array([[f2m(Rr[i, j]) for j in range(N)] for i in range(rk)], dtype=object)
    e = np.array([f2m(Rr[i, N]) for i in range(rk)], dtype=object)
    assert all(any(E[i, j] != 0 for j in range(N)) for i in range(rk)), "inconsistent equalities"
    log(f"    equalities: {rk} independent")
    # exact zeros
    rng = np.random.default_rng(31)
    xs = []
    for kind in ['k1', 'k1', 'k3', 'k3', 'rot', 'rot', 'rot', 'rot']:
        u, v = zero_vectors(kind, rng)
        assert sum(u * v) == 0 and sum(u * u) == sum(v * v), "bad zero vectors"
        L = {'u': u, 'v': v}; hats = {}
        for c in FAMS:
            m, npr = moments_exact([L[ch] for ch in c], 3)
            hats[c] = symmetrise_sites(rel, fourier_from_moments(rel, m))
        x, ok = rel.x_from_hats(hats)
        mres, ires, obj, mins = rel.evaluate(x)
        assert ok and obj == 0 and all(r == 0 for r in mres) and all(r == 0 for r in ires), f"zero {kind} failed (obj {obj})"
        xs.append(x)
        log(f"    exact zero {kind}: objective 0, moments/identities exact, min block eig {min(mins.values()):.1e}; {time.time() - t0:.0f}s")
    xstar = sum(xs[1:], xs[0])
    blocks = []; dims = [0, 0]
    for blk in rel.psd:
        name, size, d = blk
        Bx = rel.block_at(blk, xstar)
        Nk = fm2np(nullspace_q(np2fm(Bx)))
        dims[0] += size; dims[1] += Nk.shape[1]
        if Nk.shape[1] == 0: continue
        # drop blocks whose compression is identically zero (vacuous constraint)
        Nf = np.array(Nk, dtype=float)
        if max(np.abs(Nf.T @ np.array(Mi, dtype=float) @ Nf).max() for Mi in d.values()) < 1e-12: continue
        blocks.append((name, Nk.shape[1], Nk, d))
    log(f"[2] {len(xs)} exact zeros; face: block dimension {dims[0]} -> {dims[1]}, {len(blocks)} blocks kept; {time.time() - t0:.0f}s")
    # float solve of the reduced problem
    Ef = np.array(E, dtype=float); ef = np.array(e, dtype=float); cf = np.array(rel.c, dtype=float)
    esc = np.linalg.norm(Ef, axis=1); Efs = Ef / esc[:, None]; efs = ef / esc        # unit equality rows
    rows = [Efs]; b = list(efs); cones = [clarabel.ZeroConeT(Ef.shape[0])]
    nn = []; ps = []; redf = []; tr = []          # tr: per block (R^{-1} of the kernel QR, block scale)
    for name, k, Nk, d in blocks:
        Nf = np.array(Nk, dtype=float); Qk, Rk = np.linalg.qr(Nf)
        Mf = np.zeros((k * (k + 1) // 2, N)); Mexact_f = np.zeros((k * (k + 1) // 2, N))
        for i, Mi in d.items():
            Mm = np.array(Mi, dtype=float)
            Mf[:, i] = svec(Qk.T @ Mm @ Qk); Mexact_f[:, i] = svec(Nf.T @ Mm @ Nf)
        s_k = np.abs(Mf).max(); Mf = Mf / s_k
        redf.append(Mexact_f); tr.append((np.linalg.inv(Rk), s_k))
        if k == 1: nn.append(-Mf[0])
        else: ps.append((k, -Mf))
    if nn: rows.append(np.array(nn)); b += [0.0] * len(nn); cones.append(clarabel.NonnegativeConeT(len(nn)))
    for k, Mf in ps: rows.append(Mf); b += [0.0] * Mf.shape[0]; cones.append(clarabel.PSDTriangleConeT(k))
    Ad = np.vstack(rows); bvec = np.array(b)
    # parameter directions untouched by every reduced constraint are removed (the objective must vanish on them)
    ev_, V_ = np.linalg.eigh(Ad.T @ Ad); null = ev_ < 1e-13 * ev_[-1]
    if np.any(null): log(f"    Gram eigenvalues around the null threshold: largest null {ev_[null].max() / ev_[-1]:.1e}, smallest kept {ev_[~null].min() / ev_[-1]:.1e}")
    cnull = np.linalg.norm(V_[:, null].T @ cf) if np.any(null) else 0.0
    log(f"    null directions of the reduced constraints: {int(null.sum())} (objective component on them {cnull:.1e})")
    # whitening of the remaining directions: x = W xt with W = V diag(ev^-1/2), so that (Ad W)^T (Ad W) = I
    W = V_[:, ~null] / np.sqrt(ev_[~null])
    Amat = sp.csc_matrix(Ad @ W)
    import os
    st = clarabel.DefaultSettings(); st.verbose = bool(os.environ.get('PILOT_VERBOSE')); st.max_iter = 400
    st.tol_gap_abs = 1e-9; st.tol_gap_rel = 1e-9; st.tol_feas = 1e-9
    st.chordal_decomposition_enable = False
    sol = clarabel.DefaultSolver(sp.csc_matrix((W.shape[1], W.shape[1])), W.T @ cf, Amat, bvec, cones, st).solve()
    z = np.array(sol.z)
    Amat = sp.csc_matrix(Ad)
    neq = Ef.shape[0]; y0 = -z[:neq] / esc; row = neq; Sf = []
    for j in range(len(nn)): Sf.append(np.array([[z[row + j]]]))
    row += len(nn)
    for k, Mf in ps:
        kk = Mf.shape[0]; Sf.append(unsvec(z[row:row + kk], k)); row += kk
    order_idx = [i for i, blk in enumerate(blocks) if blk[1] == 1] + [i for i, blk in enumerate(blocks) if blk[1] > 1]
    order = [blocks[i] for i in order_idx]; redf_order = [redf[i] for i in order_idx]
    mineig = min(np.linalg.eigvalsh(S).min() for S in Sf)
    # back to the exact kernel coordinates: N S_N N^T = (1/s) Q S Q^T, Q = N R^{-1}  =>  S_N = (1/s) R^{-1} S R^{-T}
    S0 = []
    for S, oi in zip(Sf, order_idx):
        Rinv, s_k = tr[oi]; S0.append(Rinv @ S @ Rinv.T / s_k)
    log(f"[3] float solve: status {sol.status}, primal {sol.obj_val:.3e}, dual {-(bvec @ z):.3e}, "
        f"||A^T z + c|| {np.abs(Amat.T @ z + cf).max():.1e}, smallest reduced-dual eigenvalue {mineig:.2e}; {time.time() - t0:.0f}s")
    if mineig <= 1e-9:
        print("n=3 t=+1: FAILED (reduced dual not strictly feasible: face not minimal)", flush=True); return
    # rounding: unknowns y (neq) and the upper triangles of the S blocks
    SCALE = 2 ** 40
    rq = lambda v: mpq(int(round(float(v) * SCALE)), SCALE)
    yq = [rq(v) for v in y0]
    Sq = [[[rq(S[a, b_]) for b_ in range(S.shape[0])] for a in range(S.shape[0])] for S in S0]
    for S in Sq:
        for a in range(len(S)):
            for b_ in range(a): S[a][b_] = S[b_][a]
    def residual():
        """exact residual of c = E^T y + sum_k M_k^*(N_k S_k N_k^T) (N entries) and of e . y = 0"""
        r = np.array(rel.c, dtype=object) - E.T.dot(np.array(yq, dtype=object))
        for (name, k, Nk, d), S in zip(order, Sq):
            Sh = Nk.dot(np.array(S, dtype=object)).dot(Nk.T)
            for i, Mi in d.items(): r[i] -= (Sh * Mi).sum()
        return r, -sum(e[i] * yq[i] for i in range(neq))
    r, r_obj = residual()
    log(f"[4] rounded dual: identity residual {max(abs(float(v)) for v in r):.1e}, dual objective residual {float(r_obj):.1e}; {time.time() - t0:.0f}s")
    # column subset: all y's plus S entries chosen by pivoted QR on the float equation matrix
    cols = []          # (kind, ...) and float column vectors of length N + 1
    Fcols = []
    for i in range(neq):
        v = np.concatenate([Ef[i], [ef[i]]]); cols.append(('y', i)); Fcols.append(v)
    for bi, ((name, k, Nk, d), Mf) in enumerate(zip(order, redf_order)):
        pairs = [(a, b_) for b_ in range(k) for a in range(b_ + 1)]
        for pi, (a, b_) in enumerate(pairs):
            w = 1.0 if a == b_ else np.sqrt(2.0)
            v = np.concatenate([Mf[pi] * w / (1.0 if a == b_ else np.sqrt(2.0)) * (1.0 if a == b_ else 2.0) / (1.0 if a == b_ else np.sqrt(2.0)), [0.0]])
            cols.append(('S', bi, a, b_)); Fcols.append(v)
    Fm = np.array(Fcols).T          # (N+1) x ncols
    Qf, Rf, piv = scipy.linalg.qr(Fm, pivoting=True, mode='economic')
    rank = int(np.sum(np.abs(np.diag(Rf)) > 1e-10 * abs(Rf[0, 0])))
    chosen = [cols[j] for j in piv[:rank]]
    log(f"    column subset: equation rank {rank} of {N + 1}; {len(chosen)} unknowns corrected")
    # exact columns of the chosen unknowns
    def exact_col(cd):
        v = zeros_o((N + 1,))
        if cd[0] == 'y':
            for j in range(N): v[j] = E[cd[1], j]
            v[N] = e[cd[1]]
        else:
            _, bi, a, b_ = cd; name, k, Nk, d = order[bi]
            for i, Mi in d.items():
                val = Nk[:, a].dot(Mi).dot(Nk[:, b_])
                v[i] = val * (1 if a == b_ else 2)
        return v
    Ccols = np.array([exact_col(cd) for cd in chosen], dtype=object).T      # (N+1) x rank
    rvec = np.concatenate([r, [r_obj]])
    # independent rows of Ccols (rank many) and exact solve
    Cq = np2fm(Ccols); RR, rk2 = Cq.transpose().rref()
    rowsel = []; rr = 0
    for j in range(N + 1):
        if rr < rk2 and RR[rr, j] == 1 and all(RR[ii, j] == 0 for ii in range(rk2) if ii != rr): rowsel.append(j); rr += 1
    Csq = np2fm(Ccols[rowsel, :]); rsq = np2fm(rvec[rowsel].reshape(-1, 1))
    delta = (Csq.transpose() * (Csq * Csq.transpose()).solve(rsq)) if Csq.nrows() < Csq.ncols() else Csq.solve(rsq)
    dl = [f2m(delta[i, 0]) for i in range(delta.nrows())]
    for cd, dv in zip(chosen, dl):
        if cd[0] == 'y': yq[cd[1]] += dv
        else:
            _, bi, a, b_ = cd; Sq[bi][a][b_] += dv
            if a != b_: Sq[bi][b_][a] += dv
    r, r_obj = residual()
    exact_ok = all(v == 0 for v in r) and r_obj == 0
    log(f"[5] exact correction: max |delta| {max(abs(float(v)) for v in dl):.1e}; identity exact {exact_ok}; {time.time() - t0:.0f}s")
    if not exact_ok:
        print(f"n=3 t=+1: FAILED (identity not exact after the correction: max residual {max(abs(float(v)) for v in r):.1e})", flush=True); return
    allpsd = True; minpiv = None
    for (name, k, Nk, d), S in zip(order, Sq):
        ok, mp = exact_psd(S)
        if not ok: allpsd = False; log(f"    block {name}: NOT PSD"); break
        if mp is not None and (minpiv is None or mp < minpiv): minpiv = mp
    nb = max(max(int(abs(v.numerator)).bit_length() for v in yq), max(int(abs(x.numerator)).bit_length() for S in Sq for rw in S for x in rw))
    db = max(max(int(v.denominator).bit_length() for v in yq), max(int(x.denominator).bit_length() for S in Sq for rw in S for x in rw))
    log(f"[6] exact LDL^T: all {len(Sq)} blocks PSD {allpsd} (smallest pivot {float(minpiv) if minpiv else 0:.2e}); "
        f"rational sizes up to {nb} / {db} bits; {time.time() - t0:.0f}s")
    # independent check on an exact entangled true point (integer vectors in Z^27, normalised exactly)
    rng2 = np.random.default_rng(97)
    p = rng2.integers(-2, 3, size=27); q = rng2.integers(-2, 3, size=27)
    u = p.astype(object); v = ((p @ p) * q - (p @ q) * p).astype(object)
    L = {'u': u, 'v': v}; hats = {}
    for c in FAMS:
        m, npr = moments_exact([L[ch] for ch in c], 3)
        hats[c] = symmetrise_sites(rel, fourier_from_moments(rel, m))
    xP, okP = rel.x_from_hats(hats)
    Pm = np.outer(u, u) * mpq(1, int(u @ u)) + np.outer(v, v) * mpq(1, int(v @ v))
    q3 = exact_qn(Pm, 3); lhs = rel.c.dot(xP)
    Ex = E.dot(xP); eqok = all(Ex[i] == e[i] for i in range(neq))
    terms = []
    for (name, k, Nk, d), S in zip(order, Sq):
        Sh = Nk.dot(np.array(S, dtype=object)).dot(Nk.T)
        Mx = rel.block_at((name, Nk.shape[0], d), xP)
        terms.append((Sh * Mx).sum())
    eye_ = sum(yq[i] * Ex[i] for i in range(neq))
    log(f"[7] check on an entangled true point: in variable space {okP}; equalities hold {eqok}; c.x(P) == q_3(P): {lhs == q3} "
        f"({float(lhs):.10f}); c.x(P) == y.E x(P) + sum_k Tr(S_k M_k(x(P))): {lhs == eye_ + sum(terms)}; y.E x(P) = {eye_}; "
        f"all {len(terms)} terms >= 0: {all(tm >= 0 for tm in terms)}; {time.time() - t0:.0f}s")
    checked = okP and eqok and lhs == q3 and lhs == eye_ + sum(terms) and all(tm >= 0 for tm in terms)
    if allpsd and exact_ok and checked:
        print(f"n=3 t=+1: exact identity c = E^T y + sum_k M_k^*(N_k S_k N_k^T), e.y = 0, {len(Sq)} PSD blocks (exact LDL^T), "
              f"max {nb}/{db} bits -> CERTIFIED", flush=True)
    else:
        print("n=3 t=+1: FAILED (a reduced dual block is not PSD after rounding)", flush=True)

if __name__ == "__main__":
    main()
