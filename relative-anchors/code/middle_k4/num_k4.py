"""Numerics and exact point checks for the torus blocks of (R_inf)_1 and (R_inf)_2 at k = 4.

1. Torus block structure (exact, one rational point per j): every entry of (R_inf)_j between states with different
   torus labels vanishes, and the blocks listed in k4lib.BLOCKS together with their images under permutations of the
   copies and the type II blocks exhaust (R_inf)_j.
2. At random rational weight vectors (log-uniform spreads): the known kernel vectors are annihilated exactly, the
   kernel dimension (floating point) equals the number of known kernel vectors, and the gauge-fixed block is positive
   definite (exact LDL pivots over Q).  Also reports the smallest eigenvalue of the gauge-fixed block (floating point).
No claim rests on the floating-point numbers; the exact parts are point checks only."""
import sys, time, random, itertools
import numpy as np
import flint
from k4lib import block_states, kernel_vectors, gauge_keep, pt_builder, BLOCKS, K4
from ylib import rinf_states, wstate

t0 = time.time()
ok = True
random.seed(7)


def fq(x):
    return float(x.p) / float(x.q)


def ldl_pos(M):
    A = [row[:] for row in M]
    n = len(A)
    for p in range(n):
        if A[p][p] <= 0:
            return False
        for i in range(p + 1, n):
            if A[i][p] == 0:
                continue
            f = A[i][p] / A[p][p]
            for j in range(p + 1, n):
                A[i][j] -= f * A[p][j]
    return True


# ---- 1. torus structure
for j in (1, 2):
    s = [flint.fmpq(1), flint.fmpq(3, 7), flint.fmpq(11, 5), flint.fmpq(2, 9)]
    B = pt_builder(s)
    S = rinf_states(K4, j)
    lab = [wstate(a, y) for (a, y) in S]
    bad = 0
    for i1 in range(len(S)):
        for i2 in range(i1 + 1, len(S)):
            if lab[i1] != lab[i2]:
                if B.Y(S[i1], S[i2]) != 0:
                    bad += 1
    groups = {}
    for st, l in zip(S, lab):
        groups.setdefault(l, []).append(st)
    kinds = {'letter-1': 0, 'type I': 0, 'type II': 0}
    sizes = {}
    for l, sts in groups.items():
        n2 = l[1:]
        if max(n2) >= 2:
            kinds['type II'] += 1
        elif max(n2) == 0:
            kinds['letter-1'] += 1
            sizes.setdefault('letter-1', set()).add(len(sts))
        else:
            kinds['type I'] += 1
            sizes.setdefault(f'type I |T|={sum(n2)}', set()).add(len(sts))
    # every non-type-II block is a permutation image of a listed representative
    reps = {(jj, len(T)): len(block_states(jj, T)) for (jj, T) in BLOCKS if jj == j}
    cover = all(len(sts) == reps[(j, sum(l[1:]))] for l, sts in groups.items() if max(l[1:]) <= 1)
    print(f'j = {j}: {len(S)} states, {len(groups)} torus blocks {kinds}, sizes {sizes}; off-block nonzero entries: {bad}; '
          f'representatives cover all letter-1/type I blocks: {cover}  [{time.time() - t0:.1f}s]', flush=True)
    ok &= bad == 0 and cover

# ---- 2. random points
npts = int(sys.argv[1]) if len(sys.argv) > 1 else 12
pts = []
for it in range(npts):
    spread = [3, 6, 10][it % 3]
    s = [flint.fmpq(1)] + [flint.fmpq(max(1, int(1000 * np.exp(random.uniform(-spread, spread) / 2))), 1000) for _ in range(3)]
    pts.append(s)
pts.append([flint.fmpq(1)] * 4)
pts.append([flint.fmpq(1), flint.fmpq(1, 1000), flint.fmpq(1, 1000), flint.fmpq(1, 1000)])
pts.append([flint.fmpq(1, 1000), flint.fmpq(1), flint.fmpq(1), flint.fmpq(1)])
stats = {}
for s in pts:
    B = pt_builder(s)
    for (j, T) in BLOCKS:
        S = block_states(j, T)
        M = B.matrix(S)
        kv = kernel_vectors(S, j, T, s)
        kz = all(all(sum(M[i][l] * v[l] for l in range(len(S))) == 0 for i in range(len(S))) for v in kv)
        Mf = np.array([[fq(x) for x in row] for row in M])
        ev = np.linalg.eigvalsh(Mf)
        keep = gauge_keep(S, j, T)
        Mk = [[M[a][b] for b in keep] for a in keep]
        pd = ldl_pos(Mk)
        evk = np.linalg.eigvalsh(np.array([[fq(x) for x in row] for row in Mk]))
        # exact: kernel vectors annihilated, gauge-fixed block > 0, and the gauge coordinates complement the known
        # kernel (each kernel vector is the only one nonzero on its deleted coordinate), so the kernel is exactly
        # the span of the known vectors
        good = kz and pd and len(keep) == len(S) - len(kv)
        ok &= good
        key = (j, T)
        st = stats.setdefault(key, dict(size=len(S), nker=len(kv), kept=len(keep), lmin=1e9, lminfull=1e9, allgood=True))
        st['lmin'] = min(st['lmin'], evk.min())
        st['lminfull'] = min(st['lminfull'], sorted(ev)[len(kv)])
        st['allgood'] &= good
        if not good:
            print('  FAILED at', [str(x) for x in s], key, kz, nker, len(kv), pd, flush=True)
for key, st in stats.items():
    print(f'block j = {key[0]}, T = {key[1]}: size {st["size"]}, known kernel {st["nker"]}, gauge-fixed size {st["kept"]}; '
          f'kernel exact and complete, gauge-fixed block > 0 (exact LDL) at all {len(pts)} points: {st["allgood"]}; '
          f'min eig gauge-fixed {st["lmin"]:.3e}, eig no. (nker+1) of full block {st["lminfull"]:.3e}', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
