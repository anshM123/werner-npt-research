"""Exact check of the closed form (djlib.block) against the exact builder ylib (Lemma B, all column strings):
  2^{k-j} Y_ylib = M  entrywise, and 2^{k-j} A_ylib = A-part, at random rational weights, for every pair (k, j)
  with 1 <= j <= 6, j + 1 <= k <= 7, and (k, j) = (8, 1), (8, 3), (8, 5), (8, 7); up to three sets T each.  Also: B -> (1 + 1/1000) B in ylib equals the
  scaled closed form (used by the negative controls)."""
import sys, random, time
import flint
from ylib import Builder, point_ring
from djlib import block, states
Q = flint.fmpq
random.seed(20261105)
t0 = time.time()
ok = True; cnt = 0
pairs = [(k, j) for k in range(2, 8) for j in range(1, 7) if j + 1 <= k] + [(8, 1), (8, 3), (8, 5), (8, 7)]
for k, j in pairs:
    s = [Q(random.randint(1, 30), random.randint(1, 30)) for _ in range(k)]
    S = sum(s, Q(0)); x = [v / S for v in s]
    Bd = Builder(point_ring(s))
    Ts = set()
    for trial in range(3 if k <= 7 else 1):
        Ts.add(tuple(sorted(random.sample(range(k), j + 1))))
    for T in sorted(Ts):
        U = [u for u in range(k) if u not in T]
        sts = states(k, list(T))
        Y = [[Bd.Y(p, q) for q in sts] for p in sts]
        Ay = [[Bd.A(p[0], p[1], q[0], q[1]) for q in sts] for p in sts]
        M = block([x[t] for t in T], [x[u] for u in U])
        unit = Q(2) ** (k - j)
        good = all(unit * Y[i][l] == M[i][l] for i in range(j + 2) for l in range(j + 2))
        Mneg = block([x[t] for t in T], [x[u] for u in U], bscale=Q(1001, 1000))
        Yneg = [[Ay[i][l] - Q(1001, 1000) * (Ay[i][l] - Y[i][l]) for l in range(j + 2)] for i in range(j + 2)]
        good2 = all(unit * Yneg[i][l] == Mneg[i][l] for i in range(j + 2) for l in range(j + 2))
        cnt += 1; ok &= good and good2
        print(f'k={k} j={j} T={list(T)}: 2^(k-j) Y = closed form: {good}; scaled B: {good2} [{time.time() - t0:.1f}s]', flush=True)
print(cnt, 'blocks'); print('PASS' if ok else 'FAIL')
