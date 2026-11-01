"""Independent exact verification of the n = 4 certificate (reads only lin4.pkl, facedata4r.pkl, cert4.pkl, pts/p_ent.pkl).

Checks (all in exact rational arithmetic):
  (a) e . y = 0 and c = E^T y + sum_k F_k^*(S_k) for every coordinate;
  (b) exact LDL^T of every S_k (all pivots >= 0, zero pivots only with zero rows);
  (c) at the entangled true point P (vectors stored in pts/p_ent.pkl): q_4(P) recomputed here from P, c . x(P) = q_4(P),
      E x(P) = e, c . x(P) = y . E x(P) + sum_k Tr(S_k F_k(x(P))) with every term >= 0, and every face block
      F_k(x(P)) PSD (exact LDL^T; a necessary condition for the constraint data to be valid at a true point).
usage: python verify4.py [point_name ...]   (default: p_ent; reads facedata4r.pkl.gz if facedata4r.pkl is absent)"""
import sys, os, time, gzip, pickle, numpy as np
from gmpy2 import mpq
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'relaxations'))
from exact_coreN import exact_qn

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
    lin = pickle.load(open('lin4.pkl', 'rb')); fd = pickle.load(open('facedata4r.pkl', 'rb') if os.path.exists('facedata4r.pkl') else gzip.open('facedata4r.pkl.gz', 'rb')); cert = pickle.load(open('cert4.pkl', 'rb'))
    N = lin['N']; c = lin['c']; E = lin['E']; e = lin['e']; blocks = fd['blocks']
    y = [mpq(a, d) for (a, d) in cert['y']]
    S = {k: [[mpq(a, d) for (a, d) in row] for row in M] for k, M in cert['S'].items()}
    print(f"loaded: {N} parameters, {len(E)} equalities, {len(blocks)} face blocks, {len(S)} dual blocks; {time.time() - t0:.0f}s", flush=True)
    ey = sum((e[j] * y[j] for j in range(len(E))), mpq(0))
    r = [mpq(v) for v in c]
    for j in range(len(E)):
        if y[j] == 0: continue
        for i, v in enumerate(E[j]):
            if v != 0: r[i] -= v * y[j]
    for k, M in S.items():
        for i, v in adjoint(blocks[k], M).items(): r[i] -= v
    ident = all(v == 0 for v in r)
    print(f"(a) e . y = 0: {ey == 0};  c = E^T y + sum_k F_k^*(S_k) exactly: {ident}; {time.time() - t0:.0f}s", flush=True)
    psd = True; minp = None; sizes = []
    for k, M in S.items():
        ok, mp = ldl(M); sizes.append(len(M))
        if not ok: psd = False; print(f"    S block {blocks[k]['name']} NOT PSD"); break
        if mp is not None and (minp is None or mp < minp): minp = mp
    nb = max(int(abs(v.numerator)).bit_length() for M in S.values() for row in M for v in row)
    db = max(int(v.denominator).bit_length() for M in S.values() for row in M for v in row)
    nby = max(int(abs(v.numerator)).bit_length() for v in y); dby = max(int(v.denominator).bit_length() for v in y)
    print(f"(b) exact LDL^T: all {len(S)} dual blocks PSD: {psd} (sizes up to {max(sizes)}, smallest pivot {float(minp):.3e}); "
          f"S entries up to {nb}/{db} bits, y up to {nby}/{dby} bits; {time.time() - t0:.0f}s", flush=True)
    good = ey == 0 and ident and psd
    for pname in (sys.argv[1:] or ['p_ent']):
        rec = pickle.load(open(f'pts/{pname}.pkl', 'rb'))
        u = np.array(rec['u'], dtype=object); v = np.array(rec['v'], dtype=object)
        assert (u * v).sum() == 0
        P = np.outer(u, u) * mpq(1, int((u * u).sum())) + np.outer(v, v) * mpq(1, int((v * v).sum()))
        q4 = exact_qn(P, 4)
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
        print(f"(c) {pname}: q_4(P) recomputed from P = {float(q4):.12f}; c.x(P) == q_4(P): {cx == q4}; E x(P) == e: {eq}; "
              f"identity at the point: {cx == yEx + sum(terms)}; all {len(terms)} terms >= 0: {all(t >= 0 for t in terms)}; "
              f"all {len(blocks)} face blocks PSD at the point (exact LDL^T): {facepsd}; {time.time() - t0:.0f}s", flush=True)
        good = good and cx == q4 and eq and cx == yEx + sum(terms) and all(t >= 0 for t in terms) and facepsd
    print("VERIFIED" if good else "VERIFICATION FAILED", flush=True)

if __name__ == "__main__":
    main()
