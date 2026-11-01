"""Certificate for q_4(P) >= 0 (rank-2 projectors on (C^3)^(x)4, t = +1) from the exact face-reduced data.

Certificate: rational y (equalities E x = e) and S_k >= 0 (face blocks F_k) with
     c = E^T y + sum_k F_k^*(S_k)   exactly,   e . y = 0   exactly,
where F_k^*(S)_i = Tr(S F_k(e_i)).  Then for every rank-2 projector P (x(P) feasible): q_4(P) = c . x(P) =
y . e + sum_k Tr(S_k F_k(x(P))) >= 0.
Steps: (1) float solve of the face-reduced problem (Clarabel, whitened variables, chordal decomposition off);
(2) rounding of y and S_k to dyadic grids; (3) exact correction on a column subset of full row rank (pivoted QR on the
float equation matrix, exact square solve with flint); (4) exact identity check; (5) exact LDL^T of every S_k;
(6) independent check at an exact entangled true point.  Saves cert4.pkl.
usage: python cert4.py"""
import sys, os, time, pickle, math, itertools, numpy as np, flint, scipy.sparse as sp, scipy.linalg, clarabel
from gmpy2 import mpq, mpz

def log(s): print(s, flush=True)

def svec_pairs(m): return [(i, j) for j in range(m) for i in range(j + 1)]
def svec(V):
    m = V.shape[-1]
    return np.stack([V[..., i, j] * (1.0 if i == j else np.sqrt(2.0)) for (i, j) in svec_pairs(m)], axis=-1)
def unsvec(v, m):
    V = np.zeros((m, m)); k = 0
    for j in range(m):
        for i in range(j + 1):
            if i == j: V[i, i] = v[k]
            else: V[i, j] = V[j, i] = v[k] / np.sqrt(2.0)
            k += 1
    return V

def exact_ldl(S):
    """exact LDL^T of a symmetric rational matrix (list of lists of mpq); returns (psd, smallest positive pivot, #zero pivots)"""
    k = len(S); M = [row[:] for row in S]; minpiv = None; nzero = 0
    for i in range(k):
        p = M[i][i]
        if p < 0: return False, p, nzero
        if p == 0:
            if any(M[i][j] != 0 for j in range(i + 1, k)): return False, p, nzero
            nzero += 1; continue
        minpiv = p if minpiv is None or p < minpiv else minpiv
        for j in range(i + 1, k):
            if M[j][i] != 0:
                f = M[j][i] / p
                for l in range(i + 1, k): M[j][l] -= f * M[i][l]
    return True, minpiv, nzero

def main():
    t0 = time.time()
    with open('lin4.pkl', 'rb') as fh: lin = pickle.load(fh)
    with open('facedata4r.pkl', 'rb') as fh: fd = pickle.load(fh)
    with open('face4.pkl', 'rb') as fh: face = pickle.load(fh)
    N = lin['N']; c = lin['c']; E = lin['E']; e = lin['e']; neq = len(E)
    blocks = fd['blocks']
    log(f"[0] data: {N} parameters, {neq} equalities, {len(blocks)} face blocks "
        f"({sum(1 for b in blocks if b['kind'] == 'COV')} COV, {sum(1 for b in blocks if b['kind'] == 'X')} X), "
        f"max size {max(b['r'] for b in blocks)}; {time.time() - t0:.0f}s")
    # ------------------------------------------------------------ float data of the blocks: list over blocks of
    # (r, list of (idx array, float array (K, r, r)))
    def block_float(b):
        out = []
        if b['kind'] == 'COV':
            for gi, idx, Fh, den in b['data']:
                out.append((np.array(idx), np.array([[[float(mpq(int(v), den)) for v in row] for row in M] for M in Fh])))
        else:
            gi, idx, mats = b['data'][0]
            out.append((np.array(idx), np.array([[[a_ / b_ for (a_, b_) in row] for row in M] for M in mats], dtype=float)))
        return out
    bf = [block_float(b) for b in blocks]
    log(f"    float block data built; {time.time() - t0:.0f}s")
    # normalisation per block: R from QR of a float basis (identity coordinates) and scale s_k
    Ef = np.array([[float(v) for v in row] for row in E]); ef = np.array([float(v) for v in e]); cf = np.array([float(v) for v in c])
    esc = np.linalg.norm(Ef, axis=1); Efs = Ef / esc[:, None]; efs = ef / esc
    rows = [Efs]; bvec = list(efs); cones = [clarabel.ZeroConeT(neq)]
    nn_rows = []; ps_rows = []; trans = []          # trans[k] = (Rinv, s_k)
    for k, (b, dat) in enumerate(zip(blocks, bf)):
        r = b['r']
        # metric: average of |F_i| magnitudes -> diagonal Jacobi scaling from sum_i F_i^2 diagonal
        dg = np.zeros(r)
        for idx, F in dat: dg += np.einsum('kab,kab->b', F, F)
        assert np.all(dg > 0), f"block {b['name']} has a vacuous direction"
        dg = np.sqrt(dg)
        Rinv = np.diag(1.0 / dg)
        Mf = np.zeros((r * (r + 1) // 2, N))
        for idx, F in dat:
            G = np.einsum('ab,kbc,cd->kad', Rinv, F, Rinv)
            Mf[:, idx] = svec(G).T
        s_k = np.abs(Mf).max(); Mf /= s_k
        trans.append((Rinv, s_k))
        if r == 1: nn_rows.append((k, -Mf[0]))
        else: ps_rows.append((k, r, -Mf))
    order = [k for k, _ in nn_rows] + [k for k, _, _ in ps_rows]
    if nn_rows:
        rows.append(np.array([v for _, v in nn_rows])); bvec += [0.0] * len(nn_rows); cones.append(clarabel.NonnegativeConeT(len(nn_rows)))
    for k, r, Mf in ps_rows:
        rows.append(Mf); bvec += [0.0] * Mf.shape[0]; cones.append(clarabel.PSDTriangleConeT(r))
    Ad = np.vstack(rows); bvec = np.array(bvec); del rows
    log(f"[1] float problem: {Ad.shape[0]} rows x {N} columns ({len(nn_rows)} scalar blocks, {len(ps_rows)} PSD blocks); {time.time() - t0:.0f}s")
    Gm = Ad.T @ Ad
    ev_, V_ = np.linalg.eigh(Gm); null = ev_ < 1e-13 * ev_[-1]
    cnull = np.linalg.norm(V_[:, null].T @ cf) if np.any(null) else 0.0
    log(f"    null directions of the reduced constraints: {int(null.sum())} (objective component on them {cnull:.1e})")
    W = V_[:, ~null] / np.sqrt(ev_[~null])
    Amat = sp.csc_matrix(Ad @ W)
    st = clarabel.DefaultSettings(); st.verbose = bool(os.environ.get('CERT_VERBOSE')); st.max_iter = 400
    st.tol_gap_abs = 1e-10; st.tol_gap_rel = 1e-10; st.tol_feas = 1e-10
    st.chordal_decomposition_enable = False
    t1 = time.time()
    sol = clarabel.DefaultSolver(sp.csc_matrix((W.shape[1], W.shape[1])), W.T @ cf, Amat, bvec, cones, st).solve()
    z = np.array(sol.z)
    y0 = -z[:neq] / esc; row = neq; Sf = {}
    for j, (k, _) in enumerate(nn_rows): Sf[k] = np.array([[z[row + j]]])
    row += len(nn_rows)
    for k, r, Mf in ps_rows:
        kk = Mf.shape[0]; Sf[k] = unsvec(z[row:row + kk], r); row += kk
    mineig = min(np.linalg.eigvalsh(S).min() for S in Sf.values())
    resid = np.abs(Ad.T @ z + cf).max()
    log(f"[2] float solve: status {sol.status}, iterations {sol.iterations}, primal {sol.obj_val:.3e}, dual {-(bvec @ z):.3e}, "
        f"||A^T z + c|| {resid:.1e}, smallest reduced-dual eigenvalue {mineig:.3e}; solve {time.time() - t1:.0f}s; {time.time() - t0:.0f}s")
    del Amat, Gm
    if mineig <= 1e-9:
        log("n=4 t=+1: FAILED (reduced dual not strictly feasible: face not minimal or solve inaccurate)"); return
    # back to the exact face coordinates: Tr(S Rinv F Rinv)/s = Tr((Rinv S Rinv / s) F)
    S0 = {k: trans[k][0] @ Sf[k] @ trans[k][0].T / trans[k][1] for k in Sf}
    with open('float_dual4.pkl', 'wb') as fh: pickle.dump(dict(y0=y0, S0=S0, mineig=mineig), fh)
    # ------------------------------------------------------------ rounding to dyadic grids (per block, relative)
    def rq(v, ex):
        """v rounded to the dyadic grid 2^-ex (ex may be negative)"""
        if ex >= 0: return mpq(int(round(float(v) * 2.0 ** ex)), 2 ** ex)
        return mpq(int(round(float(v) / 2.0 ** (-ex))) * 2 ** (-ex))
    assert np.all(np.isfinite(y0)) and all(np.all(np.isfinite(S)) for S in S0.values()), "non-finite float dual"
    yex = 40 - int(np.floor(np.log2(np.abs(y0).max() + 1e-300)))
    yq = [rq(v, yex) for v in y0]
    Sq = {}
    for k, S in S0.items():
        ex = 40 - int(np.floor(np.log2(np.abs(S).max() + 1e-300)))
        r = S.shape[0]; M = [[None] * r for _ in range(r)]
        for a_ in range(r):
            for b_ in range(a_, r):
                M[a_][b_] = M[b_][a_] = rq(S[a_, b_], ex)
        Sq[k] = M
    # exact adjoint: F_k^*(S)_i for all i (object arithmetic, integer numerators)
    def adjoint(k, S):
        out = {}
        b = blocks[k]; r = b['r']
        den_S = 1
        for row in S:
            for v in row: den_S = den_S * int(v.denominator) // math.gcd(den_S, int(v.denominator))
        St = np.array([[int(v * den_S) for v in row] for row in S], dtype=object).ravel()
        if b['kind'] == 'COV':
            for gi, idx, Fh, den in b['data']:
                vals = Fh.reshape(Fh.shape[0], r * r).dot(St)
                for kk, i in enumerate(idx): out[i] = out.get(i, mpq(0)) + mpq(int(vals[kk]), den * den_S)
        else:
            gi, idx, mats = b['data'][0]
            for kk, i in enumerate(idx):
                tot = mpq(0)
                for a_ in range(r):
                    for b_ in range(r):
                        num, dd = mats[kk][a_][b_]
                        if num: tot += mpq(num, dd) * S[a_][b_]
                out[i] = out.get(i, mpq(0)) + tot
        return out
    def residual():
        rres = [mpq(v) for v in c]
        for j in range(neq):
            yj = yq[j]
            if yj == 0: continue
            for i, v in enumerate(E[j]):
                if v != 0: rres[i] -= v * yj
        for k in Sq:
            for i, v in adjoint(k, Sq[k]).items(): rres[i] -= v
        robj = -sum((e[j] * yq[j] for j in range(neq)), mpq(0))
        return rres, robj
    rres, robj = residual()
    log(f"[3] rounded dual: identity residual {max(abs(float(v)) for v in rres):.1e}, dual objective residual {float(robj):.1e}; {time.time() - t0:.0f}s")
    # ------------------------------------------------------------ column subset
    cols = []; Fc = []
    for j in range(neq):
        v = np.concatenate([Ef[j], [ef[j]]]); cols.append(('y', j)); Fc.append(v / np.linalg.norm(v))
    for k in Sq:
        r = blocks[k]['r']; dat = bf[k]
        for (a_, b_) in svec_pairs(r):
            v = np.zeros(N + 1)
            for idx, F in dat: v[idx] += F[:, a_, b_] * (1.0 if a_ == b_ else 2.0)
            nv = np.linalg.norm(v)
            if nv == 0: continue
            cols.append(('S', k, a_, b_)); Fc.append(v / nv)
    Fm = np.array(Fc).T; del Fc
    log(f"    equation matrix {Fm.shape[0]} x {Fm.shape[1]}; pivoted QR ...")
    Qf, Rf, piv = scipy.linalg.qr(Fm, pivoting=True, mode='economic')
    dR = np.abs(np.diag(Rf)); rank = int(np.sum(dR > 1e-10 * dR[0]))
    chosen = [cols[j] for j in piv[:rank]]
    log(f"    column subset: equation rank {rank} of {N + 1} (diag R ratio at cut {dR[rank - 1] / dR[0]:.1e}, next "
        f"{(dR[rank] / dR[0]) if rank < len(dR) else 0:.1e}); {sum(1 for cd in chosen if cd[0] == 'y')} y's, "
        f"{sum(1 for cd in chosen if cd[0] == 'S')} dual entries; {time.time() - t0:.0f}s")
    del Fm, Qf, Rf
    # exact columns
    def exact_col(cd):
        v = [mpq(0)] * (N + 1)
        if cd[0] == 'y':
            for i in range(N): v[i] = E[cd[1]][i]
            v[N] = e[cd[1]]
        else:
            _, k, a_, b_ = cd; b = blocks[k]; r = b['r']; f = 1 if a_ == b_ else 2
            if b['kind'] == 'COV':
                for gi, idx, Fh, den in b['data']:
                    for kk, i in enumerate(idx):
                        val = Fh[kk, a_, b_]
                        if val: v[i] += mpq(int(val) * f, den)
            else:
                gi, idx, mats = b['data'][0]
                for kk, i in enumerate(idx):
                    num, dd = mats[kk][a_][b_]
                    if num: v[i] += mpq(num * f, dd)
        return v
    C = [exact_col(cd) for cd in chosen]          # list of columns
    rvec = rres + [robj]
    # integer scaling of columns for flint; row selection by float QR of the chosen columns
    Cf = np.array([[float(x) for x in col] for col in C])          # rank x (N+1)
    _, _, rpiv = scipy.linalg.qr(Cf, pivoting=True, mode='economic')
    rowsel = sorted(int(j) for j in rpiv[:rank])
    cden = []
    for col in C:
        d = 1
        for x in col: d = d * int(x.denominator) // math.gcd(d, int(x.denominator))
        cden.append(d)
    Asq = flint.fmpz_mat(rank, rank, [int(C[jj][ii] * cden[jj]) for ii in rowsel for jj in range(rank)])
    bsq = flint.fmpq_mat(rank, 1, [flint.fmpq(int(rvec[ii].numerator), int(rvec[ii].denominator)) for ii in rowsel])
    t1 = time.time()
    Aq = flint.fmpq_mat(Asq)
    sol_d = Aq.solve(bsq, algorithm='dixon')
    dl = [mpq(int(sol_d[jj, 0].p), int(sol_d[jj, 0].q)) * cden[jj] for jj in range(rank)]
    log(f"    exact square solve ({rank} x {rank}) in {time.time() - t1:.0f}s; correction bits up to "
        f"{max(max(int(abs(v.numerator)).bit_length(), int(v.denominator).bit_length()) for v in dl)}")
    for cd, dv in zip(chosen, dl):
        if cd[0] == 'y': yq[cd[1]] += dv
        else:
            _, k, a_, b_ = cd; Sq[k][a_][b_] += dv
            if a_ != b_: Sq[k][b_][a_] += dv
    rres, robj = residual()
    exact_ok = all(v == 0 for v in rres) and robj == 0
    log(f"[4] exact correction: max |delta| {max(abs(float(v)) for v in dl):.1e}; identity c = E^T y + sum_k F_k^*(S_k) exact: "
        f"{exact_ok}; e . y = 0: {robj == 0}; {time.time() - t0:.0f}s")
    with open('cert4.pkl', 'wb') as fh:
        pickle.dump(dict(y=[(int(v.numerator), int(v.denominator)) for v in yq],
                         S={k: [[(int(v.numerator), int(v.denominator)) for v in row] for row in M] for k, M in Sq.items()}), fh)
    if not exact_ok:
        log(f"n=4 t=+1: FAILED (identity not exact after the correction: max residual {max(abs(float(v)) for v in rres):.1e})"); return
    # ------------------------------------------------------------ exact LDL^T
    allpsd = True; minpiv = None; t1 = time.time()
    for k, M in Sq.items():
        ok, mp, nz = exact_ldl(M)
        if not ok: allpsd = False; log(f"    block {blocks[k]['name']}: NOT PSD"); break
        if mp is not None and (minpiv is None or mp < minpiv): minpiv = mp
    nb = max(max(int(abs(v.numerator)).bit_length() for v in yq), max(int(abs(x.numerator)).bit_length() for M in Sq.values() for rw in M for x in rw))
    db = max(max(int(v.denominator).bit_length() for v in yq), max(int(x.denominator).bit_length() for M in Sq.values() for rw in M for x in rw))
    log(f"[5] exact LDL^T: all {len(Sq)} blocks PSD {allpsd} (smallest pivot {float(minpiv) if minpiv else 0:.2e}); rational sizes up to "
        f"{nb} / {db} bits; LDL {time.time() - t1:.0f}s; {time.time() - t0:.0f}s")
    # ------------------------------------------------------------ independent check at the entangled true point
    with open('pts/p_ent.pkl', 'rb') as fh: rec = pickle.load(fh)
    xP = [mpq(p_, q_) for (p_, q_) in rec['x']]
    q4 = mpq(*rec['q4'])
    lhs = sum((c[i] * xP[i] for i in range(N) if c[i] != 0), mpq(0))
    Ex = [sum((E[j][i] * xP[i] for i in range(N) if E[j][i] != 0), mpq(0)) for j in range(neq)]
    eqok = all(Ex[j] == e[j] for j in range(neq))
    yEx = sum((yq[j] * Ex[j] for j in range(neq)), mpq(0))
    terms = []
    for k, M in Sq.items():
        tot = mpq(0)
        for i, v in adjoint(k, M).items():
            if xP[i] != 0: tot += v * xP[i]
        terms.append(tot)
    checked = rec['in_space'] and eqok and lhs == q4 and lhs == yEx + sum(terms) and all(tm >= 0 for tm in terms)
    log(f"[6] entangled true point: in variable space {rec['in_space']}; equalities hold {eqok}; c.x(P) == q_4(P): {lhs == q4} "
        f"({float(lhs):.10f}); c.x(P) == y.E x(P) + sum_k Tr(S_k F_k(x(P))): {lhs == yEx + sum(terms)}; y.E x(P) = {yEx}; "
        f"all {len(terms)} terms >= 0: {all(tm >= 0 for tm in terms)} (smallest {float(min(terms)):.2e}); {time.time() - t0:.0f}s")
    if allpsd and exact_ok and checked:
        log(f"n=4 t=+1: exact identity c = E^T y + sum_k F_k^*(S_k), e.y = 0, {len(Sq)} PSD blocks (exact LDL^T), "
            f"max {nb}/{db} bits, entangled true-point check passed -> CERTIFIED")
    else:
        log("n=4 t=+1: FAILED (" + ("block not PSD" if not allpsd else "true-point check failed") + ")")

if __name__ == "__main__":
    main()
