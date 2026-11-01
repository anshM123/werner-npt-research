"""Audit copy of verifyL.py: more ratios (including near-degenerate ones) and the near-zero audit points by default.

Independent exact verification of the lifted n = 4 certificate: q_4(a P_u + b P_v) >= 0 for all real a, b.

Reads only linL.pkl, facedataLr.pkl, the certificate (default certL.pkl) and the point files ../exact4/pts/*.pkl,
../exact4m/pts/*.pkl; imports no pipeline code (numpy, gmpy2 and the standard library only).
Checks (exact rational arithmetic throughout):
  (0) the saved objectives satisfy c+ = (fu + fv + 2 fuv)/16, c- = (fu + fv - 2 fuv)/16, d = (fu - fv)/16;
  (a) e . y+ = e . y- = e . yZ = 0 and, in every coordinate i,
          c+_i/4 = (E^T y+)_i + sum_k Tr(A_k F+_k(e_i)),   c-_i/4 = (E^T y-)_i + sum_k Tr(C_k F-_k(e_i)),
          d_i/2  = (E^T yZ)_i + 2 sum_k Tr(B_k^T G_k(e_i)),          S_k = [[A_k, B_k], [B_k^T, C_k]];
  (b) exact LDL^T of every S_k (pivots >= 0, a zero pivot only with a zero row);
  (c) at each entangled true point (u, v real integer vectors, x = x(u, v) from the point file):
      u . v = 0, E x = e, every joint block H_k(x) = [[F+_k(x), G_k(x)], [G_k(x)^T, F-_k(x)]] PSD (exact fraction-free
      elimination of the integer matrix (common denominator) H_k(x); then
      L_k = diag((a+b) I, (a-b) I) H_k(x) diag((a+b) I, (a-b) I) is PSD for every ratio), and for each ratio (a, b):
      q_4(a P_u + b P_v) recomputed here from the 81 x 81 matrix (own integer partial traces), the lifted objective
      (1/4)[(a+b)^2 c+ . x + (a-b)^2 c- . x + 2 (a^2 - b^2) d . x] equal to it, and the identity
          objective = sum_g w_g y_g . E x + sum_k Tr(S_k L_k)   (w = (a+b)^2, (a-b)^2, a^2 - b^2)
      with y_g . E x = 0 and every term Tr(S_k L_k) >= 0.
usage: python verifyL.py [--cert FILE] [point ...]   (point = exact4/<name> or exact4m/<name>; default: the entangled
       points p_ent, p_ent2, a_ent1..4 of ../exact4/pts and p_ent, p_ent2, b_ent1..3 of ../exact4m/pts)"""
import sys, os, time, gzip, pickle, itertools, math, numpy as np
from gmpy2 import mpq

RATIOS = [(1, 0), (0, 1), (1, 1), (1, -1), (1000, 999), (1000, -999), (7, -2), (13, 29), (-17, 4), (2, 1), (1, -3)]
DEFAULT_BUILD = ['exact4/p_ent', 'exact4/p_ent2', 'exact4/a_ent1', 'exact4/a_ent2', 'exact4/a_ent3', 'exact4/a_ent4',
           'exact4m/p_ent', 'exact4m/p_ent2', 'exact4m/b_ent1', 'exact4m/b_ent2', 'exact4m/b_ent3']
DEFAULT = ['exact4/a_near1', 'exact4/a_near4', 'exact4m/b_nsf1', 'exact4m/b_nsf2', 'exact4m/b_nk21', 'exact4m/b_nin1',
           'exact4/a_ent1', 'exact4m/b_ent2']

def ldl(S):
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

def psd_int(M):
    """exact PSD test of a symmetric integer matrix (list of lists of Python ints) by fraction-free symmetric elimination
    (Bareiss): after the positive pivots P the active entries are the minors det M[P+j, P+l] = det M[P, P] * (Schur
    complement)_jl, so signs and zeros are those of the rational LDL^T; a zero pivot must come with a zero row (then the
    index is skipped); every division is exact (Sylvester's identity) and is asserted"""
    n = len(M); A = [list(row) for row in M]; prev = 1
    for i in range(n):
        p = A[i][i]
        if p < 0: return False
        if p == 0:
            if any(A[i][j] != 0 for j in range(i + 1, n)): return False
            continue
        Ai = A[i]
        for j in range(i + 1, n):
            Aj = A[j]; aji = Aj[i]
            for l in range(j, n):
                q, rem = divmod(p * Aj[l] - aji * Ai[l], prev)
                assert rem == 0, "inexact Bareiss division"
                Aj[l] = q
            for l in range(j + 1, n): A[l][j] = Aj[l]
        prev = p
    return True

def q4_ratio(u, v, a, b):
    """exact q_4(a P_u + b P_v) = sum_J (-1/2)^|J| ||Tr_J (a P_u + b P_v)||^2 for real integer u, v (81 entries) and integers
    a, b, via the integer matrix M = |u|^2 |v|^2 (a P_u + b P_v) = a |v|^2 u u^T + b |u|^2 v v^T: for each J the traced sites
    are moved to the end and the partial trace is the sum of the diagonal slices"""
    nu = int(sum(t * t for t in u)); nv = int(sum(t * t for t in v)); D = nu * nv
    M = np.empty((81, 81), dtype=object)
    for i in range(81):
        for j in range(81): M[i, j] = a * nv * int(u[i]) * int(u[j]) + b * nu * int(v[i]) * int(v[j])
    tot = mpq(0)
    for k in range(5):
        for J in itertools.combinations(range(4), k):
            keep = [s for s in range(4) if s not in J]
            T = M.reshape((3,) * 8).transpose(keep + [4 + s for s in keep] + list(J) + [4 + s for s in J])
            T = T.reshape(3 ** len(keep), 3 ** len(keep), 3 ** k, 3 ** k)
            R = T[:, :, 0, 0].copy()
            for i in range(1, 3 ** k): R = R + T[:, :, i, i]
            tot += mpq(-1, 2) ** k * sum(int(z) * int(z) for z in R.ravel())
    return tot / (D * D)

def joint_at(b, aint):
    """integer matrix L * H_k(x) (x = aint / lambda, L * lambda > 0 a common scale) and the scale L (the lcm of the D's)"""
    rp, rm = b['rp'], b['rm']; r = rp + rm
    L = 1
    for (_, _, _, _, _, D) in b['data']: L = L * D // math.gcd(L, D)
    H = np.zeros((r, r), dtype=object)
    for gi, idx, Fp, Fm, G, D in b['data']:
        av = np.array([aint[i] for i in idx], dtype=object)
        if not any(t != 0 for t in av): continue
        f = L // D
        if rp: H[:rp, :rp] += np.tensordot(av, Fp, axes=1) * f
        if rm: H[rp:, rp:] += np.tensordot(av, Fm, axes=1) * f
        if rp and rm:
            Gx = np.tensordot(av, G, axes=1) * f
            H[:rp, rp:] += Gx; H[rp:, :rp] += Gx.T
    return H, L

def main():
    t0 = time.time()
    args = sys.argv[1:]; certfile = 'certL.pkl'
    if args and args[0] == '--cert': certfile = args[1]; args = args[2:]
    points = args or DEFAULT
    lin = pickle.load(open('linL.pkl', 'rb')); fd = pickle.load(open('facedataLr.pkl', 'rb') if os.path.exists('facedataLr.pkl') else gzip.open('facedataLr.pkl.gz', 'rb')); cert = pickle.load(open(certfile, 'rb'))
    N = lin['N']; E = lin['E']; e = lin['e']; neq = len(E); blocks = fd['blocks']
    cp, cm, dd, fu, fv, fuv = lin['cplus'], lin['cminus'], lin['d'], lin['fu'], lin['fv'], lin['fuv']
    ok0 = all(cp[i] == (fu[i] + fv[i] + 2 * fuv[i]) / 16 and cm[i] == (fu[i] + fv[i] - 2 * fuv[i]) / 16 and
              dd[i] == (fu[i] - fv[i]) / 16 for i in range(N))
    y = {g: [mpq(a, q) for (a, q) in cert['y'][g]] for g in ('+', '-', 'Z')}
    S = {k: [[mpq(a, q) for (a, q) in row] for row in M] for k, M in cert['S'].items()}
    assert set(S) == set(range(len(blocks))), "dual blocks do not match the data blocks"
    for k, M in S.items():
        r = blocks[k]['rp'] + blocks[k]['rm']
        assert len(M) == r and all(len(row) == r for row in M) and all(M[i][j] == M[j][i] for i in range(r) for j in range(r))
    print(f"loaded: {N} parameters, {neq} equalities, {len(blocks)} lifted blocks (total dimension "
          f"{sum(b['rp'] + b['rm'] for b in blocks)}, largest {max(b['rp'] + b['rm'] for b in blocks)}), certificate {certfile}; "
          f"(0) c+, c-, d consistent with fu, fv, fuv: {ok0}; {time.time() - t0:.0f}s", flush=True)
    # ---------------------------------------------------------------- (a) identities
    tgt = {'+': [v / 4 for v in cp], '-': [v / 4 for v in cm], 'Z': [v / 2 for v in dd]}
    good = ok0
    for g in ('+', '-', 'Z'):
        res = list(tgt[g])
        for j in range(neq):
            if y[g][j] == 0: continue
            for i, v in enumerate(E[j]):
                if v != 0: res[i] -= v * y[g][j]
        for k, b in enumerate(blocks):
            rp, rm = b['rp'], b['rm']; Sa = np.array(S[k], dtype=object)
            part = {'+': Sa[:rp, :rp], '-': Sa[rp:, rp:], 'Z': Sa[:rp, rp:]}[g]
            if part.size == 0: continue
            for gi, idx, Fp, Fm, G, D in b['data']:
                Md = {'+': Fp, '-': Fm, 'Z': G}[g]
                if Md.size == 0: continue
                vals = Md.reshape(Md.shape[0], -1).dot(part.ravel())
                f = 2 if g == 'Z' else 1
                for kk, i in enumerate(idx): res[i] -= vals[kk] * f / D
        ey = sum((e[j] * y[g][j] for j in range(neq)), mpq(0))
        nbad = sum(1 for v in res if v != 0)
        good = good and ey == 0 and nbad == 0
        print(f"(a) identity {g}: e . y{g} = 0: {ey == 0}; exact in all {N} coordinates: {nbad == 0} ({nbad} nonzero residuals); "
              f"{time.time() - t0:.0f}s", flush=True)
    # ---------------------------------------------------------------- (b) PSD
    psd = True; minp = None
    for k, M in S.items():
        ok, mp = ldl(M)
        if not ok: psd = False; print(f"    S block {k} ({blocks[k]['name']}) NOT PSD", flush=True); break
        if mp is not None and (minp is None or mp < minp): minp = mp
    nb = max(int(abs(v.numerator)).bit_length() for M in S.values() for row in M for v in row)
    db = max(int(v.denominator).bit_length() for M in S.values() for row in M for v in row)
    nby = max(int(abs(v.numerator)).bit_length() for g in y for v in y[g]); dby = max(int(v.denominator).bit_length() for g in y for v in y[g])
    good = good and psd
    print(f"(b) exact LDL^T: all {len(S)} lifted dual blocks PSD: {psd} (smallest pivot {float(minp) if minp is not None else 0:.3e}); "
          f"S entries up to {nb}/{db} bits, y up to {nby}/{dby} bits; {time.time() - t0:.0f}s", flush=True)
    # ---------------------------------------------------------------- (c) true points
    npass = 0; ncases = 0
    for pname in points:
        folder, nm = pname.split('/')
        rec = pickle.load(open(f'../{folder}/pts/{nm}.pkl', 'rb'))
        u = [int(t) for t in rec['u']]; v = [int(t) for t in rec['v']]
        orth = sum(s * t for s, t in zip(u, v)) == 0 and any(u) and any(v)
        x = [mpq(a, q) for (a, q) in rec['x']]
        lam = 1
        for t in x: lam = lam * int(t.denominator) // math.gcd(lam, int(t.denominator))
        aint = [int(t * lam) for t in x]
        Ex = [sum((E[j][i] * x[i] for i in range(N) if E[j][i] != 0), mpq(0)) for j in range(neq)]
        eq = all(Ex[j] == e[j] for j in range(neq))
        yEx = {g: sum((y[g][j] * Ex[j] for j in range(neq)), mpq(0)) for g in y}
        dot = lambda c: sum((c[i] * x[i] for i in range(N) if c[i] != 0), mpq(0))
        cpx, cmx, dx = dot(cp), dot(cm), dot(dd)
        tA = []; tB = []; tC = []; jpsd = True
        for k, b in enumerate(blocks):
            H, L = joint_at(b, aint)
            r = b['rp'] + b['rm']; rp = b['rp']
            ok = psd_int([[int(H[i, j]) for j in range(r)] for i in range(r)])
            if not ok: jpsd = False; print(f"    joint block {b['name']} NOT PSD at {pname}", flush=True)
            sc = mpq(1, L * lam)
            Sk = S[k]
            tA.append(sum((Sk[i][j] * int(H[i, j]) for i in range(rp) for j in range(rp) if H[i, j] != 0), mpq(0)) * sc)
            tC.append(sum((Sk[i][j] * int(H[i, j]) for i in range(rp, r) for j in range(rp, r) if H[i, j] != 0), mpq(0)) * sc)
            tB.append(sum((Sk[i][j] * int(H[i, j]) for i in range(rp) for j in range(rp, r) if H[i, j] != 0), mpq(0)) * sc)
        ptok = orth and eq and jpsd
        lines = []
        for (a, b_) in RATIOS:
            wp, wm, wz = (a + b_) ** 2, (a - b_) ** 2, a * a - b_ * b_
            q4 = q4_ratio(u, v, a, b_)
            obj = (wp * cpx + wm * cmx + 2 * wz * dx) / 4
            terms = [wp * tA[k] + 2 * wz * tB[k] + wm * tC[k] for k in range(len(blocks))]
            yt = wp * yEx['+'] + wm * yEx['-'] + wz * yEx['Z']
            ok = obj == q4 and obj == yt + sum(terms) and yt == 0 and all(t >= 0 for t in terms)
            ncases += 1; npass += ok; ptok = ptok and ok
            lines.append(f"({a},{b_}): q4 {float(q4):.6e} obj==q4 {obj == q4} id {obj == yt + sum(terms)} y {yt == 0} "
                         f"terms>=0 {all(t >= 0 for t in terms)}")
        good = good and ptok
        print(f"(c) {pname}: u.v = 0 {orth}; E x = e {eq}; all {len(blocks)} joint blocks H_k(x) PSD {jpsd}; "
              f"{time.time() - t0:.0f}s\n      " + "\n      ".join(lines), flush=True)
    print(f"(c) summary: {npass} of {ncases} (point, ratio) cases passed", flush=True)
    print("VERIFIED" if good else "VERIFICATION FAILED", flush=True)

if __name__ == "__main__":
    main()
