"""Independent exact verification of the n = 4, t = -1 certificate (reads only lin4m.pkl, facedata4mr.pkl, the
certificate file (default cert4m.pkl) and pts/<point>.pkl; no pipeline code is imported).

Statement checked: q_4(P_u - P_v) >= 0 on the relaxation, i.e. c = E^T y + sum_k F_k^*(S_k) with e . y = 0 and S_k >= 0,
c the t = -1 objective.  Checks (all in exact rational arithmetic):
  (a) e . y = 0 and c = E^T y + sum_k F_k^*(S_k) for every coordinate;
  (b) every S_k is keyed to an existing face block of the right size, is exactly symmetric, and passes an exact LDL^T
      (all pivots >= 0, zero pivots only with zero rows);
  (c) at each entangled true point (vectors u, v stored in pts/<name>.pkl): q_4(P_u - P_v) recomputed here from the
      81 x 81 matrix P = P_u - P_v (own partial-trace code, integer arithmetic), c . x(P) = q_4(P_u - P_v), E x(P) = e,
      c . x(P) = y . E x(P) + sum_k Tr(S_k F_k(x(P))) with every term >= 0, and every face block F_k(x(P)) PSD (exact
      LDL^T; a necessary condition for the constraint data to be valid at a true point).
usage: python verify4m.py [--cert FILE] [point_name ...]   (default point: p_ent)"""
import sys, os, time, gzip, pickle, itertools, numpy as np
from gmpy2 import mpq

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

def q4_minus(u, v):
    """exact q_4(P_u - P_v) = sum_J (-1/2)^|J| ||Tr_J P||^2 for real integer u, v (81 entries), via the integer matrix
    M = |u|^2 |v|^2 P = |v|^2 u u^T - |u|^2 v v^T: for each J the traced sites are moved to the end and the partial
    trace is the sum of the diagonal (i, i) slices"""
    nu = int(sum(a * a for a in u)); nv = int(sum(a * a for a in v)); D = nu * nv
    M = np.empty((81, 81), dtype=object)
    for i in range(81):
        for j in range(81): M[i, j] = nv * int(u[i]) * int(u[j]) - nu * int(v[i]) * int(v[j])
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

def block_at(b, x):
    """exact F_k(x) as list of lists of mpq"""
    r = b['r']; tot = np.full((r, r), mpq(0), dtype=object)
    if b['kind'] == 'COV':
        for gi, idx, Fh, den in b['data']:
            xv = np.array([x[i] for i in idx], dtype=object)
            tot = tot + np.tensordot(xv, Fh, axes=1) * mpq(1, den)
    else:
        gi, idx, mats = b['data'][0]
        for k, i in enumerate(idx):
            if x[i] != 0: tot = tot + np.array([[mpq(a, d) for (a, d) in row] for row in mats[k]], dtype=object) * x[i]
    return [[tot[a, b_] for b_ in range(r)] for a in range(r)]

def adjoint(b, S):
    r = b['r']; out = {}
    Sa = np.array(S, dtype=object)
    if b['kind'] == 'COV':
        for gi, idx, Fh, den in b['data']:
            vals = Fh.reshape(Fh.shape[0], r * r).dot(Sa.ravel())
            for kk, i in enumerate(idx): out[i] = out.get(i, mpq(0)) + vals[kk] * mpq(1, den)
    else:
        gi, idx, mats = b['data'][0]
        for kk, i in enumerate(idx):
            M = np.array([[mpq(a, d) for (a, d) in row] for row in mats[kk]], dtype=object)
            out[i] = out.get(i, mpq(0)) + (M * Sa).sum()
    return out

def main():
    t0 = time.time()
    args = sys.argv[1:]; certfile = 'cert4m.pkl'
    if args and args[0] == '--cert': certfile = args[1]; args = args[2:]
    points = args or ['p_ent']
    lin = pickle.load(open('lin4m.pkl', 'rb')); fd = pickle.load(open('facedata4mr.pkl', 'rb') if os.path.exists('facedata4mr.pkl') else gzip.open('facedata4mr.pkl.gz', 'rb')); cert = pickle.load(open(certfile, 'rb'))
    assert lin['t'] == -1
    N = lin['N']; c = lin['c']; E = lin['E']; e = lin['e']; blocks = fd['blocks']
    y = [mpq(a, d) for (a, d) in cert['y']]
    S = {k: [[mpq(a, d) for (a, d) in row] for row in M] for k, M in cert['S'].items()}
    print(f"loaded: {N} parameters, {len(E)} equalities, {len(blocks)} face blocks, {len(S)} dual blocks (certificate "
          f"{certfile}); {time.time() - t0:.0f}s", flush=True)
    ey = sum((e[j] * y[j] for j in range(len(E))), mpq(0))
    r = [mpq(v) for v in c]
    for j in range(len(E)):
        if y[j] == 0: continue
        for i, v in enumerate(E[j]):
            if v != 0: r[i] -= v * y[j]
    for k, M in S.items():
        for i, v in adjoint(blocks[k], M).items(): r[i] -= v
    ident = all(v == 0 for v in r)
    nbad = sum(1 for v in r if v != 0)
    print(f"(a) e . y = 0: {ey == 0};  c = E^T y + sum_k F_k^*(S_k) exactly: {ident} ({nbad} nonzero residual coordinates); "
          f"{time.time() - t0:.0f}s", flush=True)
    psd = True; minp = None; sizes = []
    keys_ok = all(isinstance(k, int) and 0 <= k < len(blocks) and len(M) == blocks[k]['r'] and all(len(row) == len(M) for row in M)
                  for k, M in S.items())
    sym = all(M[a][b_] == M[b_][a] for M in S.values() for a in range(len(M)) for b_ in range(a + 1, len(M)))
    for k, M in S.items():
        ok, mp = ldl(M); sizes.append(len(M))
        if not ok: psd = False; print(f"    S block {blocks[k]['name']} NOT PSD"); break
        if mp is not None and (minp is None or mp < minp): minp = mp
    nb = max(int(abs(v.numerator)).bit_length() for M in S.values() for row in M for v in row)
    db = max(int(v.denominator).bit_length() for M in S.values() for row in M for v in row)
    nby = max(int(abs(v.numerator)).bit_length() for v in y); dby = max(int(v.denominator).bit_length() for v in y)
    print(f"(b) dual blocks keyed to face blocks of the right size: {keys_ok}; exactly symmetric: {sym}; "
          f"exact LDL^T: all {len(S)} dual blocks PSD: {psd} (sizes up to {max(sizes)}, smallest pivot "
          f"{float(minp) if minp is not None else 0:.3e}); S entries up to {nb}/{db} bits, y up to {nby}/{dby} bits; "
          f"{time.time() - t0:.0f}s", flush=True)
    good = ey == 0 and ident and psd and keys_ok and sym
    for pname in points:
        rec = pickle.load(open(f'pts/{pname}.pkl', 'rb'))
        u = [int(a) for a in rec['u']]; v = [int(a) for a in rec['v']]
        orth = sum(a * b_ for a, b_ in zip(u, v)) == 0 and any(u) and any(v)
        q4 = q4_minus(u, v)
        x = [mpq(a, d) for (a, d) in rec['x']]
        cx = sum((c[i] * x[i] for i in range(N) if c[i] != 0), mpq(0))
        Ex = [sum((E[j][i] * x[i] for i in range(N) if E[j][i] != 0), mpq(0)) for j in range(len(E))]
        eq = all(Ex[j] == e[j] for j in range(len(E)))
        yEx = sum((y[j] * Ex[j] for j in range(len(E))), mpq(0))
        terms = []; facepsd = True
        for k, b in enumerate(blocks):
            Fx = block_at(b, x)
            ok, _ = ldl(Fx)
            if not ok: facepsd = False; print(f"    face block {b['name']} NOT PSD at the true point")
            if k in S: terms.append(sum((S[k][a][b_] * Fx[a][b_] for a in range(b['r']) for b_ in range(b['r'])), mpq(0)))
        print(f"(c) {pname}: u.v = 0: {orth}; q_4(P_u - P_v) recomputed from P = {float(q4):.12f} (stored value agrees: "
              f"{mpq(*rec['q4']) == q4}); c.x(P) == q_4(P_u - P_v): {cx == q4}; E x(P) == e: {eq}; identity at the point: "
              f"{cx == yEx + sum(terms)}; all {len(terms)} terms >= 0: {all(t >= 0 for t in terms)}; all {len(blocks)} face "
              f"blocks PSD at the point (exact LDL^T): {facepsd}; {time.time() - t0:.0f}s", flush=True)
        good = good and orth and cx == q4 and eq and cx == yEx + sum(terms) and all(t >= 0 for t in terms) and facepsd
    print("VERIFIED" if good else "VERIFICATION FAILED", flush=True)

if __name__ == "__main__":
    main()
