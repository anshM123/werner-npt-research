"""Cross-check of the modular compression + CRT of facedata4.py against a direct exact recomputation.

For a random sample of covariance face blocks and every parameter group touching them, recompute
    Fhat = F0 + F0^T,  F0 = sum_{r,s} (Gt U N_r)^T Ytil_rs (U N_s)
with Python integers only (no primes), from the float64-exact tree images, and compare with facedata4.pkl entry by
entry; also check that facedata4r.pkl holds exactly B Fhat B^T for the same blocks.
usage: python crosscheck4.py [nblocks]"""
import sys, time, pickle, numpy as np
from core4 import *
from engine4 import *

def main():
    nsample = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    t0 = time.time(); n = 4
    A, B1, A3 = site_algebras()
    vs = VarSpace(A, n, log=lambda s: None)
    face = pickle.load(open('face4.pkl', 'rb')); fd = pickle.load(open('facedata4.pkl', 'rb')); fdr = pickle.load(open('facedata4r.pkl', 'rb'))
    mus = sorted_patterns(len(B1.blocks), n)
    Gts = {mu: np.array(fm_to_int(fm_kron([B1.blocks[q]['G'] for q in mu]), 1), dtype=object) for mu in mus}
    info = {}
    for (xy, mu), cl in face['U'].items():
        for lab, dsig, U in cl:
            Nk = face['cker'][(xy, mu, lab)]; m = U.shape[1]
            if Nk.shape[1] == 0: continue
            NL = [Nk[:m], Nk[m:]]
            info[f"COV{xy}{mu}{lab}"] = (xy, mu, [Gts[mu].dot(U).dot(NL[r]) for r in range(2)], [U.dot(NL[s_]) for s_ in range(2)])
    maps = {}
    for (xy, r, s_, c, P) in cov_entries():
        Pm = PERM[P]
        maps[(xy, r, s_)] = (c,) + site_map_int(A, B1, lambda X, Pm=Pm: realign1(Pm @ X @ Pm.T))
    tau = max(rowsum_bound(m[1]) for m in maps.values())
    cov = [b for b in fd['blocks'] if b['kind'] == 'COV']
    rng = np.random.default_rng(2026)
    pick = sorted(rng.choice(len(cov), size=min(nsample, len(cov)), replace=False))
    red = {b['name']: b for b in fdr['blocks']}
    nent = 0; bad = 0; badr = 0; trees = {}
    for j in pick:
        b = cov[j]; xy, mu, L, R = info[b['name']]
        for gi, idx, Fh, den in b['data']:
            key = vs.order[gi]; g = vs.groups[key]
            F0 = np.zeros((g['K'], b['r'], b['r']), dtype=object)
            for r in range(2):
                for s_ in range(2):
                    c, Ti, dT = maps[(xy, r, s_)]
                    if c != key[0]: continue
                    tk = (xy, r, s_, key)
                    if tk not in trees: trees[tk] = tree(vs.pair_float(key), key[1], Ti, len(B1.blocks))
                    Y = assemble(trees[tk], g['orbit'], mu)
                    if Y is None: continue
                    check_exact(Y, g['Xmax'] * tau ** n * len(g['orbit']))
                    Yi = pair_to_mat(Y, tuple(B1.blocks[q]['f'] for q in mu)).astype(np.int64).astype(object)
                    for k in range(g['K']): F0[k] += L[r].T.dot(Yi[k]).dot(R[s_])
            Fd = F0 + F0.transpose(0, 2, 1)
            nent += Fd.size
            if not all(a_ == b_ for a_, b_ in zip(Fd.ravel(), Fh.ravel())): bad += 1
            # reduced block: B Fhat B^T
            rb = red[b['name']]; Bm = rb['B']
            got = [d for d in rb['data'] if d[0] == gi][0][2]
            want = np.array([Bm.dot(Fd[k]).dot(Bm.T) for k in range(g['K'])], dtype=object)
            if not all(a_ == b_ for a_, b_ in zip(want.ravel(), got.ravel())): badr += 1
        trees = {k: v for k, v in trees.items() if k[0] == xy}
    print(f"crosscheck: {len(pick)} covariance face blocks, {nent} numerator entries recomputed with Python integers: "
          f"{'all equal' if bad == 0 else f'{bad} (block, group) MISMATCHES'} to the CRT values; reduced data B Fhat B^T: "
          f"{'all equal' if badr == 0 else f'{badr} MISMATCHES'}; {time.time() - t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
