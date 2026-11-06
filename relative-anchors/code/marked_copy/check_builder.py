"""Exact checks of the reduced forms against the Lemma B builder (j1lib.py), and end-to-end kernels at small k.
(a) For k = 2..KB, every j = -1..k, every i, several rational squares sigma and bscale in {1, 1001/1000}: every entry
    of the reduced form of mlib (numeric mode) equals v_a^T Y~ v_b, where Y~ is the rescaled j1lib block and v_a the
    explicit twisted vectors of explicit.py; the norms agree.
(b) For k = 2..KE and several sigma: the known kernel vectors of the whole block (column-group indicators; the K2
    vector at j = k - 2; the k - 1 Theorem T vectors at j = k - 1) are annihilated by the j1lib block, its exact
    nullity equals their number, and it is positive semidefinite (exact symmetric elimination).
(c) The symbolic reduced forms evaluated at rational sigma equal the numeric ones (k up to KS).
(d) The Newton window: the numeric forms with all Johnson nodes equal those with the window i-2..i+2 (k up to KS).
Usage: python check_builder.py KB KE KS"""
import sys, time
from fractions import Fraction as Fr
import flint
import mlib, explicit

KB, KE, KS = (int(x) for x in sys.argv[1:4])
fq = flint.fmpq
T0 = time.time()
ok = True


def F(x):
    return Fr(int(x.p), int(x.q))


def ldl_psd(Mx):
    """exact symmetric elimination with diagonal pivoting: (is_psd, rank)"""
    n = len(Mx)
    M = [[fq(x.numerator, x.denominator) for x in row] for row in Mx]
    act = list(range(n))
    rank = 0
    while act:
        cand = [t for t in act if M[t][t] != 0]
        if not cand:
            return all(M[a][b] == 0 for a in act for b in act), rank
        t = cand[0]
        piv = M[t][t]
        if piv < 0:
            return False, rank
        rank += 1
        act.remove(t)
        for a in act:
            f = M[a][t] / piv
            if f != 0:
                for b in act:
                    M[a][b] -= f * M[t][b]
    return True, rank


# (a)
SQ = [Fr(2), Fr(1, 3), Fr(3, 2), Fr(1, 10), Fr(7, 4), Fr(5)]
cnt = 0
for k in range(2, KB + 1):
    ws = SQ if k <= 5 else SQ[:3]
    for w0 in ws:
        sig = w0 * w0
        w = [w0] + [Fr(1)] * (k - 1)
        for bs in (Fr(1), Fr(1001, 1000)):
            for j in range(-1, k + 1):
                st, M = explicit.block(k, j, w, bscale=bs)
                for i in range(0, (k - 1) // 2 + 1):
                    names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=sig, bscale=bs))
                    st2, vecs = explicit.twisted(k, j, i)
                    if sorted(names) != sorted(vecs):
                        ok = False
                        print('NAMES DIFFER', k, j, i, names, sorted(vecs))
                        continue
                    for a, na in enumerate(names):
                        if sum(c * c for c in vecs[na].values()) != norms[a]:
                            ok = False
                            print('NORM', k, j, i, na)
                        for b, nb in enumerate(names):
                            cnt += 1
                            if explicit.form(M, vecs[na], vecs[nb]) != F(Q[a][b]):
                                ok = False
                                print('ENTRY DIFFERS', k, j, i, sig, bs, na, nb)
    print(f'(a) k={k}: done ({time.time() - T0:.0f}s)', flush=True)
print(f'(a) {cnt} entries compared with the Lemma B builder', flush=True)

# (b)
npts = 0
for k in range(2, KE + 1):
    for w0 in (Fr(2), Fr(1, 3), Fr(1, 10)):
        sig = w0 * w0
        w = [w0] + [Fr(1)] * (k - 1)
        for j in range(-1, k + 1):
            st, M = explicit.block(k, j, w)
            n = len(st)
            kv = explicit.known_kernel(k, j, sig)
            for v in kv:
                for a in range(n):
                    if sum(M[a][x] * c for x, c in v.items()) != 0:
                        ok = False
                        print('KNOWN VECTOR NOT IN KERNEL', k, j, sig)
                        break
            psd, rank = ldl_psd(M)
            npts += 1
            if not psd or n - rank != len(kv):
                ok = False
                print('END TO END FAILS', k, j, sig, 'psd', psd, 'nullity', n - rank, 'expected', len(kv))
    print(f'(b) k={k}: done ({time.time() - T0:.0f}s)', flush=True)
print(f'(b) {npts} blocks: known vectors annihilated, PSD, nullity = number of known vectors', flush=True)

# (c), (d)
cnt_c = cnt_d = 0
for k in range(2, KS + 1):
    for j in range(-1, k + 1):
        for i in range(0, (k - 1) // 2 + 1):
            names, norms, Q = mlib.reduced(mlib.Calc(k, j, i, sigma=None))
            for sig in (Fr(2, 9), Fr(11, 3)):
                n2, nr2, Q2 = mlib.reduced(mlib.Calc(k, j, i, sigma=sig))
                n3, nr3, Q3 = mlib.reduced(mlib.Calc(k, j, i, sigma=sig, full_nodes=True))
                if n2 != names or n3 != names:
                    ok = False
                    print('NAMES (c/d)', k, j, i)
                    continue
                for a in range(len(names)):
                    for b in range(len(names)):
                        cnt_c += 1
                        cnt_d += 1
                        if Q[a][b](sig) != F(Q2[a][b]):
                            ok = False
                            print('SYMBOLIC != NUMERIC', k, j, i, sig)
                        if Q3[a][b] != Q2[a][b]:
                            ok = False
                            print('WINDOW != ALL NODES', k, j, i, sig)
print(f'(c) {cnt_c} symbolic entries evaluated = numeric; (d) {cnt_d} window = all nodes', flush=True)
print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')
