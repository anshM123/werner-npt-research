"""Build every torus block of (R_inf)_j, j = -1..k, for k = 2 and 3 (one representative per orbit of the copy
permutations), exactly over Q(s_0, ..., s_{k-1}) in the rescaled basis (ylib.Builder), with the greedy gauge of
poslib.gauge.  Checks (exact):
  - the symbolic block is symmetric and equals the point builder at 3 random rational points;
  - the kernel dimension is the same at 3 random rational points, and the computed kernel at one point is annihilated
    by the block at the others only through the gauge statement below;
  - the gauge-fixed principal submatrix is nonsingular (indeed positive definite, exact LDL) at 3 random points, and
    deleted + kept = all states with |deleted| = dim kernel.
Stores blocks_k{k}.pkl: list of dicts (label, kind, states, keep, dele, M = full block as (num dict, den dict)).
Usage: python build_blocks.py k"""
import sys, time, pickle, random
import flint
from poslib import representatives, label_state, sym_builder, pt_builder, gauge, rand_point, nullspace_Q, Q

k = int(sys.argv[1])
t0 = time.time()
rng = random.Random(1000 + k)
Bs, K, s = sym_builder(k, dehom=False)
pts = [rand_point(k, rng) for _ in range(3)]
Bp = [pt_builder(p) for p in pts]
out = []
ok = True


def evalF(f, p):
    return Q(f.n(*p)) / Q(f.d(*p))


def is_pd(Mx):
    n = len(Mx)
    A = [row[:] for row in Mx]
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            for l in range(c, n):
                A[i][l] -= f * A[c][l]
    return True


for j in range(-1, k + 1):
    reps, B = representatives(k, j, kinds_wanted=('L1', 'I', 'II'))
    for lab, kind in reps:
        S = B[lab]
        n = len(S)
        M = Bs.matrix(S)
        sym = all(M[a][b] == M[b][a] for a in range(n) for b in range(n))
        agree = True
        kdims = []
        gpd = []
        Mps = []
        for p, bp in zip(pts, Bp):
            Mp = bp.matrix(S)
            Mps.append(Mp)
            agree &= all(evalF(M[a][b], p) == Mp[a][b] for a in range(n) for b in range(n))
            kdims.append(len(nullspace_Q(Mp)))
        ker, dele, keep = gauge(S, Mps[0])
        for Mp in Mps:
            gpd.append(is_pd([[Mp[a][b] for b in keep] for a in keep]))
        good = sym and agree and len(set(kdims)) == 1 and all(gpd) and len(keep) + len(dele) == n
        ok &= good
        print(f'k = {k}, j = {j}, label {lab} ({kind}), size {n}, kernel dim {kdims}, gauge-fixed size {len(keep)}, '
              f'deleted {[label_state(S[i]) for i in dele]}; symmetric {sym}, agrees with point builder {agree}, '
              f'gauge-fixed PD at 3 points {gpd}  [{time.time() - t0:.1f}s]', flush=True)
        out.append(dict(k=k, j=j, label=lab, kind=kind, states=S, keep=keep, dele=dele,
                        labels=[label_state(x) for x in S],
                        M=[[(M[a][b].n.to_dict(), M[a][b].d.to_dict()) for b in range(n)] for a in range(n)]))
pickle.dump(dict(k=k, names=[f's{t}' for t in range(k)], blocks=out), open(f'blocks_k{k}.pkl', 'wb'))
print(f'{len(out)} blocks; time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
