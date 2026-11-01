"""Exact lifted span reduction (as ../exact4/span4.py, now for the lifted blocks).

Lifted block k: L_k = [[F+_k(Y+), G_k(Z)], [G_k(Z)^T, F-_k(Y-)]] on R^{d+} (+) R^{d-}.  Its data matrices are, for every
parameter i, [[F+_i, 0], [0, 0]] (Y+_i), [[0, 0], [0, F-_i]] (Y-_i) and [[0, G_i], [G_i^T, 0]] (Z_i), so all data
vanish outside V+ (+) V-,
    V+ = sum_i range(F+_i) + sum_i colspace(G_i)  = row space of the stack [F+_i ; G_i^T]  (rows in R^{d+}),
    V- = sum_i range(F-_i) + sum_i colspace(G_i^T) = row space of the stack [F-_i ; G_i]  (rows in R^{d-}).
B+ (r+ x d+), B- (r- x d-): LLL-reduced integer bases of these row spaces (rows chosen independent modulo
p = 1048573, exact RREF); completeness checked EXACTLY (every stacked matrix times an exact basis of ker B+- is 0).
Reduced data: F+'_i = B+ F+_i B+^T, F-'_i = B- F-_i B-^T, G'_i = B+ G_i B-^T, i.e. the compression
diag(B+, B-) L_k diag(B+, B-)^T, which is PSD whenever L_k is (valid for any B), and exactly equivalent here.
Inputs (read only): ../exact4/facedata4.pkl (F+), ../exact4m/facedata4m.pkl (F-), facedataL.pkl (G, d+, d-).
Saves facedataLr.pkl: blocks = list of dict(name, kind, rp, rm, Bp, Bm, data = [(gi, idx, Fp, Fm, G, D)]) with integer
object arrays Fp (K, rp, rp), Fm (K, rm, rm), G (K, rp, rm) and a common denominator D per (block, group):
F+'_k(e_i) = Fp[kk] / D, F-'_k(e_i) = Fm[kk] / D, G'_k(e_i) = G[kk] / D for i = idx[kk].
usage: python spanL.py"""
import sys, time, pickle, math, numpy as np, flint
from build_linL import independent_rows_modp, P_SEL

def frac_mats_to_int(mats):
    d = 1
    for M in mats:
        for row in M:
            for (a, q) in row: d = d * q // math.gcd(d, q)
    return np.array([[[a * (d // q) for (a, q) in row] for row in M] for M in mats], dtype=object), d

def load_face(fn):
    fd = pickle.load(open(fn, 'rb')); out = {}
    for b in fd['blocks']:
        if b['kind'] == 'COV': out[b['name']] = {gi: (idx, Fh, den) for gi, idx, Fh, den in b['data']}
        else:
            gi, idx, mats = b['data'][0]; arr, d = frac_mats_to_int(mats); out[b['name']] = {gi: (idx, arr, d)}
    return out

def row_basis(St, ncol):
    """LLL-reduced integer basis (object array rk x ncol) of the row space of the integer matrix St (object), with
    the exact completeness check St . Kb = 0 for an exact basis Kb of ker B"""
    if St.shape[0] == 0 or not any(v != 0 for v in St.ravel()): return np.zeros((0, ncol), dtype=object)
    Stp = np.array([[int(v) % P_SEL for v in row] for row in St], dtype=np.int64)
    sel = independent_rows_modp(Stp)
    rows = St[sel]
    Rq = flint.fmpq_mat(len(sel), ncol, [flint.fmpq(int(v)) for v in rows.ravel()])
    RR, rk = Rq.rref()
    assert rk == len(sel)
    d = 1
    for i in range(rk):
        for j in range(ncol): d = d * int(RR[i, j].q) // math.gcd(d, int(RR[i, j].q))
    Z = flint.fmpz_mat(rk, ncol, [int((RR[i, j] * d).p) for i in range(rk) for j in range(ncol)]).lll()
    B = np.array([[int(Z[i, j]) for j in range(ncol)] for i in range(rk)], dtype=object).reshape(rk, ncol)
    Kb, nul = flint.fmpz_mat(rk, ncol, [int(v) for v in B.ravel()]).nullspace()
    if nul:
        Kbo = np.array([[int(Kb[i, j]) for j in range(nul)] for i in range(ncol)], dtype=object).reshape(ncol, nul)
        prod = St.dot(Kbo)
        assert all(v == 0 for v in prod.ravel()), "row space incomplete modulo p"
    return B

def main():
    t0 = time.time()
    Fp_all = load_face('../exact4/facedata4.pkl'); Fm_all = load_face('../exact4m/facedata4m.pkl')
    fl = pickle.load(open('facedataL.pkl', 'rb'))
    print(f"[spanL] loaded: {len(Fp_all)} t=+1 face blocks, {len(Fm_all)} t=-1 face blocks, {len(fl['blocks'])} lifted blocks; "
          f"{time.time() - t0:.0f}s", flush=True)
    names = set(b['name'] for b in fl['blocks'])
    assert set(Fp_all) <= names and set(Fm_all) <= names
    newblocks = []; tot = [0, 0]; tot_sep = [0, 0]; nG = 0
    for b in fl['blocks']:
        name, dp, dm = b['name'], b['dp'], b['dm']
        Gd = {gi: (idx, Gh, den) for gi, idx, Gh, den in b['data']}
        Pd = Fp_all.get(name, {}); Md = Fm_all.get(name, {})
        for gi, (idx, Fh, den) in Pd.items(): assert Fh.shape[1:] == (dp, dp), name
        for gi, (idx, Fh, den) in Md.items(): assert Fh.shape[1:] == (dm, dm), name
        # stacks
        sp_ = [Fh[k] for gi, (idx, Fh, den) in Pd.items() for k in range(Fh.shape[0])] + \
              [Gh[k].T for gi, (idx, Gh, den) in Gd.items() for k in range(Gh.shape[0])]
        sm_ = [Fh[k] for gi, (idx, Fh, den) in Md.items() for k in range(Fh.shape[0])] + \
              [Gh[k] for gi, (idx, Gh, den) in Gd.items() for k in range(Gh.shape[0])]
        Bp = row_basis(np.vstack(sp_) if sp_ and dp else np.zeros((0, dp), dtype=object), dp)
        Bm = row_basis(np.vstack(sm_) if sm_ and dm else np.zeros((0, dm), dtype=object), dm)
        rp, rm = Bp.shape[0], Bm.shape[0]
        # separate spans (F+ alone, F- alone) for comparison with span4 / span4m
        sep_p = row_basis(np.vstack([Fh[k] for gi, (idx, Fh, den) in Pd.items() for k in range(Fh.shape[0])]), dp).shape[0] if Pd else 0
        sep_m = row_basis(np.vstack([Fh[k] for gi, (idx, Fh, den) in Md.items() for k in range(Fh.shape[0])]), dm).shape[0] if Md else 0
        tot[0] += dp + dm; tot[1] += rp + rm; tot_sep[0] += sep_p; tot_sep[1] += sep_m
        if Gd: nG += 1
        if rp + rm == 0: continue
        data = []
        for gi in sorted(set(Pd) | set(Md) | set(Gd)):
            idx = (Pd.get(gi) or Md.get(gi) or Gd.get(gi))[0]
            for dct in (Pd, Md, Gd):
                if gi in dct: assert list(dct[gi][0]) == list(idx)
            dens = [dct[gi][2] for dct in (Pd, Md, Gd) if gi in dct]
            D = 1
            for dd in dens: D = D * dd // math.gcd(D, dd)
            K = len(idx)
            Fp = np.zeros((K, rp, rp), dtype=object); Fm = np.zeros((K, rm, rm), dtype=object); G = np.zeros((K, rp, rm), dtype=object)
            if gi in Pd and rp:
                _, Fh, den = Pd[gi]; Fp = np.array([Bp.dot(Fh[k]).dot(Bp.T) * (D // den) for k in range(K)], dtype=object).reshape(K, rp, rp)
            if gi in Md and rm:
                _, Fh, den = Md[gi]; Fm = np.array([Bm.dot(Fh[k]).dot(Bm.T) * (D // den) for k in range(K)], dtype=object).reshape(K, rm, rm)
            if gi in Gd and rp and rm:
                _, Gh, den = Gd[gi]; G = np.array([Bp.dot(Gh[k]).dot(Bm.T) * (D // den) for k in range(K)], dtype=object).reshape(K, rp, rm)
            if not (any(v != 0 for v in Fp.ravel()) or any(v != 0 for v in Fm.ravel()) or any(v != 0 for v in G.ravel())): continue
            data.append((gi, list(idx), Fp, Fm, G, D))
        newblocks.append(dict(name=name, kind=b['kind'], rp=rp, rm=rm, Bp=Bp, Bm=Bm, data=data))
    bbits = max(max((abs(int(v)) for v in np.concatenate([b['Bp'].ravel(), b['Bm'].ravel()])), default=0).bit_length() for b in newblocks)
    nbits = max(max((abs(int(v)).bit_length() for (_, _, Fp, Fm, G, D) in b['data'] for M in (Fp, Fm, G) for v in M.ravel()), default=0)
                for b in newblocks)
    sizes = [b['rp'] + b['rm'] for b in newblocks]
    nsvec = sum(r * (r + 1) // 2 for r in sizes)
    pickle.dump(dict(blocks=newblocks), open('facedataLr.pkl', 'wb'))
    print(f"[spanL] {len(fl['blocks'])} lifted blocks ({nG} with nonzero cross data); total dimension {tot[0]} -> {tot[1]} "
          f"(separate spans of F+ / F- alone: {tot_sep[0]} / {tot_sep[1]}); {len(newblocks)} blocks kept (largest {max(sizes)}, "
          f"svec {nsvec}, {sum(1 for b in newblocks if b['rp'] and b['rm'])} with both parts); B entries <= {bbits} bits, "
          f"numerators <= {nbits} bits; row-space completeness checked exactly; {time.time() - t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
