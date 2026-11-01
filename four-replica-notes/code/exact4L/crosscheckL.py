"""Cross-check of the modular compression + CRT of facedataL.py (cross data) against a direct exact recomputation.

For a random sample of COV blocks with nonzero cross data and every parameter group touching them, recompute
    Ghat = sum_{r,s} (Gt U N+_r)^T Ytil_rs (U N-_s) + [ sum_{r,s} (Gt U N-_r)^T Ytil_rs (U N+_s) ]^T
with Python integers only (no primes), from the float64-exact tree images, and compare with facedataL.pkl entry by
entry; also check that facedataLr.pkl holds exactly B+ Ghat B-^T (times D / den) for the same blocks.
usage: python crosscheckL.py [nblocks]"""
import sys, time, pickle, numpy as np
from core4 import *
from engine4 import *

def main():
    nsample = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    t0 = time.time(); n = 4
    A, B1, A3 = site_algebras()
    vs = VarSpace(A, n, log=lambda s: None)
    fp = pickle.load(open('../exact4/face4.pkl', 'rb')); fm = pickle.load(open('../exact4m/face4m.pkl', 'rb'))
    fl = pickle.load(open('facedataL.pkl', 'rb')); flr = pickle.load(open('facedataLr.pkl', 'rb'))
    mus = sorted_patterns(len(B1.blocks), n)
    Gts = {mu: np.array(fm_to_int(fm_kron([B1.blocks[q]['G'] for q in mu]), 1), dtype=object) for mu in mus}
    info = {}
    for (xy, mu), cl in fp['U'].items():
        for lab, dsig, U in cl:
            Np = fp['cker'][(xy, mu, lab)]; Nm = fm['cker'][(xy, mu, lab)]; m = U.shape[1]
            if Np.shape[1] == 0 or Nm.shape[1] == 0: continue
            GU = Gts[mu].dot(U)
            info[f"COV{xy}{mu}{lab}"] = (xy, mu, [GU.dot(Np[r * m:(r + 1) * m]) for r in range(2)], [U.dot(Np[r * m:(r + 1) * m]) for r in range(2)],
                                         [GU.dot(Nm[r * m:(r + 1) * m]) for r in range(2)], [U.dot(Nm[r * m:(r + 1) * m]) for r in range(2)])
    maps = {}
    for (xy, r, s_, c, P) in cov_entries():
        Pm = PERM[P]
        maps[(xy, r, s_)] = (c,) + site_map_int(A, B1, lambda X, Pm=Pm: realign1(Pm @ X @ Pm.T))
    tau = max(rowsum_bound(m[1]) for m in maps.values())
    cov = [b for b in fl['blocks'] if b['kind'] == 'COV' and b['data']]
    rng = np.random.default_rng(2027)
    pick = sorted(rng.choice(len(cov), size=min(nsample, len(cov)), replace=False))
    red = {b['name']: b for b in flr['blocks']}
    nent = 0; bad = 0; badr = 0; nred = 0; trees = {}
    for j in pick:
        b = cov[j]; xy, mu, Lp, Rp, Lm, Rm = info[b['name']]
        for gi, idx, Gh, den in b['data']:
            key = vs.order[gi]; g = vs.groups[key]
            G0 = np.zeros((g['K'], b['dp'], b['dm']), dtype=object); G1 = np.zeros((g['K'], b['dm'], b['dp']), dtype=object)
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
                    for k in range(g['K']):
                        G0[k] += Lp[r].T.dot(Yi[k]).dot(Rm[s_])
                        G1[k] += Lm[r].T.dot(Yi[k]).dot(Rp[s_])
            Gd = G0 + G1.transpose(0, 2, 1)
            assert den == 2 * 8 ** n * g['dX']
            nent += Gd.size
            if not all(a_ == b_ for a_, b_ in zip(Gd.ravel(), Gh.ravel())): bad += 1
            rb = red.get(b['name'])
            if rb is None: continue
            got = [d for d in rb['data'] if d[0] == gi]
            if not got: badr += 1; continue
            _, _, Fp, Fm, Gr, D = got[0]
            want = np.array([rb['Bp'].dot(Gd[k]).dot(rb['Bm'].T) * (D // den) for k in range(g['K'])], dtype=object).reshape(Gr.shape)
            nred += 1
            if not all(a_ == b_ for a_, b_ in zip(want.ravel(), Gr.ravel())): badr += 1
        trees = {k: v for k, v in trees.items() if k[0] == xy}
    print(f"crosscheckL: {len(pick)} COV blocks with cross data, {nent} cross numerator entries recomputed with Python integers: "
          f"{'all equal' if bad == 0 else f'{bad} (block, group) MISMATCHES'} to the CRT values; reduced cross data B+ Ghat B-^T "
          f"({nred} block-groups): {'all equal' if badr == 0 else f'{badr} MISMATCHES'}; {time.time() - t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
