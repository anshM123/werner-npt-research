"""Bookkeeping of every block type behind Corollary I2.

python blocks_t.py KE KMIN KMAX PKL [PKL ...]
(a) For k = 3..KE and every j = -1..k: the states of (R_inf)_j (crosscheck.rinf_states: all rows and columns with any
    letters) are grouped by the torus label (number of letters 2 at each copy, row and column together; the number of
    letters 1 is then fixed).  Every group is the letter-1 block, a type II block (some copy with two letters 2), or a
    type I block (labels 0/1, T = the copies with label 1); for every type I block with a letter 1 the group equals
    explicit_t.block_states(k, j, T) (the state list used by the end-to-end checks).
(b) For every k in [KMIN, KMAX] and every j = -1..k, the classes of blocks and the result that covers each at squared
    weights (sigma, 1, ..., 1) are listed; for 3 <= j <= k - 1 every type I block with a letter 1 is reduced (Section 1
    of NOTES.md) to (j, |T|, case), and for each such triple every form of tlib.all_forms must have a certificate in
    the given files (status 'PD on complement')."""
import sys, pickle
from math import comb
from crosscheck import rinf_states, label
from explicit_t import block_states
import tlib


def main():
    ke, kmin, kmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    saved = {}
    for p in sys.argv[4:]:
        with open(p, 'rb') as fh:
            d = pickle.load(fh)
        for key, info in d.items():      # keep only what the bookkeeping needs (memory)
            saved[key] = {f: info.get(f) for f in ('status', 'names', 'nkernel', 'gauge', 'ngauged')}
        del d
    ok = True
    for k in range(3, ke + 1):
        cnt = {'letter-1': 0, 'type II': 0, 'type I no letter 1': 0, 'type I with a letter 1': 0}
        for j in range(-1, k + 1):
            groups = {}
            for s in rinf_states(k, j):
                groups.setdefault(label(s, k), []).append(s)
            for lab, gs in groups.items():
                if any(x == 2 for x in lab):
                    cnt['type II'] += 1
                elif all(x == 0 for x in lab):
                    cnt['letter-1'] += 1
                else:
                    T = [t for t in range(k) if lab[t] == 1]
                    if len(T) > j:
                        cnt['type I no letter 1'] += 1
                    else:
                        cnt['type I with a letter 1'] += 1
                        if set(gs) != set(block_states(k, j, T)):
                            ok = False
                            print('BLOCK STATES DIFFER', k, j, T, flush=True)
        print(f'(a) k={k}: blocks by class {cnt}; every type I block with a letter 1 equals block_states', flush=True)
    for k in range(kmin, kmax + 1):
        missing = []
        ntrip = 0
        for j in range(3, k):
            byj = {}
            for f in tlib.all_forms(k, j, j):
                byj.setdefault((f[1], f[2]), []).append(f)
            for tau in range(1, j + 1):
                for case in ('MT', 'MU'):
                    if case == 'MU' and tau > k - 1:
                        continue
                    ntrip += 1
                    forms = byj.get((tau, case), [])
                    if not forms:
                        missing.append((j, tau, case, 'no forms'))
                    for f in forms:
                        info = saved.get((k,) + f)
                        if info is None or info.get('status') != 'PD on complement':
                            missing.append(f)
        if missing:
            ok = False
            print(f'(b) k={k}: MISSING {missing[:5]} ({len(missing)})', flush=True)
        nt = sum(comb(k, t) for j in range(3, k) for t in range(1, j + 1))
        print(f'(b) k={k}: (R_inf)_-1, _0, _k: Proposition M, Theorems B0, K; (R_inf)_1: Section 3h corollary (two '
              f'weight values); (R_inf)_2: Corollary M2; for 3 <= j <= k-1: letter-1 blocks Theorem MK, blocks without '
              f'a letter 1 Theorem DJ, type II Proposition T2prime, type I with a letter 1 ({nt} labels T, reduced to '
              f'{ntrip} triples (j, |T|, case)) Theorem I1: all certified {not missing}', flush=True)
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
