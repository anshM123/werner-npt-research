"""Coverage bookkeeping of the fixed-k certificates.

python coverage_t.py KMIN KMAX PKL [PKL ...]
For every k in [KMIN, KMAX]: every (j, tau, case, mode, i) with 3 <= j <= k - 1, 1 <= tau <= j (the type I blocks
with a letter 1), the marked copy in T or in U, T'-mode sym or anti and irrep index i = 0..floor(u'/2), is enumerated
here directly; its twisted space is nonzero iff some twisted vector has a nonzero norm (tlib.norm_of, which
check_explicit.py compares with the explicit vectors for k <= 8).  Each nonzero one must have a saved certificate
with status 'PD on complement', the gauge of one coordinate per known kernel vector, and the expected number of known
kernel vectors (sym mode: one per value of the column letter at the marked copy for which the column group is
nonempty when the marked copy is in U, one when it is in T; anti mode: none).  Also printed per k: the number of type
I blocks (j, tau, case) and of forms, and that the anti mode occurs exactly when tau' >= 2."""
import sys, pickle
import tlib


def expected_kernel(k, j, tau, case, mode, i, names):
    if mode == 'anti':
        return 0
    if case == 'MT':
        return 1 if any(nm in names for nm in ('A0', 'A1', 'B10', 'B01', 'C', 'CM')) else 0
    cnt = 0
    for y0 in (0, 1):
        if any(nm in names for nm in ('O%d' % y0, 'A0%d' % y0, 'A1%d' % y0, 'B10%d' % y0, 'B01%d' % y0, 'C%d' % y0)):
            cnt += 1
    return cnt


def main():
    kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
    saved = {}
    for p in sys.argv[3:]:
        with open(p, 'rb') as fh:
            d = pickle.load(fh)
        for key, info in d.items():      # keep only what the bookkeeping needs (memory)
            saved[key] = {f: info.get(f) for f in ('status', 'names', 'nkernel', 'gauge', 'ngauged')}
        del d
    ok = True
    total = 0
    for k in range(kmin, kmax + 1):
        nblocks = 0
        nforms = {('MT', 'sym'): 0, ('MT', 'anti'): 0, ('MU', 'sym'): 0, ('MU', 'anti'): 0}
        missing = []
        for j in range(3, k):
            for tau in range(1, j + 1):
                for case in ('MT', 'MU'):
                    up = k - tau if case == 'MT' else k - tau - 1
                    if up < 0:
                        continue
                    nblocks += 1
                    for mode in ('sym', 'anti'):
                        taup = tau - 1 if case == 'MT' else tau
                        for i in range(0, up // 2 + 1):
                            C = tlib.Calc(k, j, tau, i, case, mode, sigma=1)
                            names = [nm for nm in tlib.NAMES[(case, mode)] if C.valid() and tlib.norm_of(C, nm) != 0]
                            if mode == 'anti' and taup < 2:
                                assert not C.valid()
                            if not names:
                                continue
                            nforms[(case, mode)] += 1
                            key = (k, j, tau, case, mode, i)
                            info = saved.get(key)
                            if info is None:
                                missing.append(key)
                                continue
                            ek = expected_kernel(k, j, tau, case, mode, i, names)
                            if (info.get('status') != 'PD on complement' or info.get('names') != names
                                    or info.get('nkernel') != ek or len(info.get('gauge', [])) != ek
                                    or info.get('ngauged') != len(names) - ek):
                                missing.append(key + ('bad record',))
        nf = sum(nforms.values())
        total += nf
        if missing:
            ok = False
            print(f'k={k}: MISSING or bad: {missing[:5]} ({len(missing)} in all)', flush=True)
        print(f'k={k}: {nblocks} type I blocks with a letter 1 (j = 3..{k - 1}, both cases), {nf} reduced forms '
              f'(MT sym {nforms[("MT", "sym")]}, MT anti {nforms[("MT", "anti")]}, MU sym {nforms[("MU", "sym")]}, '
              f'MU anti {nforms[("MU", "anti")]}), all certified: {not missing}', flush=True)
    print(f'{total} forms in all')
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
