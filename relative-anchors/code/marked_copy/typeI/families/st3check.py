"""A check of the family certificates for the type I blocks with a letter 1 at one marked copy (Theorems AT, S1, SD,
S2), independent of them: at sampled members (k, j, tau, case, mode, i) with 30 <= k <= KMAX that lie in a certified
family, the numeric reduced form of tlib.py at a rational sigma is positive semidefinite with kernel exactly the span of
its known kernel vectors: the kernel vectors are annihilated exactly, the rank is n - (number of kernel vectors), and
a principal submatrix of that size with nonzero determinant is positive definite (exact symmetric elimination).
python st3check.py NPER SEED KMAX PKL [PKL ...]      (NPER members per (case, mode, family) sampled at most)"""
import sys, random, pickle, time, itertools
from fractions import Fraction as Fr
import flint
import tlib

nper, seed, kmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
fams = []
for p in sys.argv[4:]:
    d = pickle.load(open(p, 'rb'))
    for key, r in d['results'].items():
        if isinstance(r, str) and r.startswith('PD'):
            fams.append(key)
rng = random.Random(seed)
sigmas = [Fr(1, 10 ** 5), Fr(1, 13), Fr(2, 7), Fr(1), Fr(9, 4), Fr(31), Fr(10 ** 5)]


def member(spec):
    """a random member (q0, t, i, tau') of a family spec, or None"""
    out = []
    for s in spec:
        if isinstance(s, tuple):
            out.append(s[1] + rng.randint(0, 30))
        else:
            out.append(s)
    return tuple(out)


def pd(M):
    n = len(M)
    A = [row[:] for row in M]
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] / A[c][c]
                for cc in range(c + 1, n):
                    A[r][cc] -= f * A[c][cc]
    return True


ok = True
T0 = time.time()
ndone = 0
for key in fams:
    case, mode = key[0], key[1]
    spec = key[2:]
    got = 0
    for attempt in range(200):
        if got >= nper:
            break
        q0, t, i, tp = member(spec)
        tau = tp + 1 if case == 'MT' else tp
        # MT: tau = tau' + 1 and k = tau + u' with u' = m + 2i; MU: tau = tau' and k = tau + u' + 1 (coverage_st3.forms)
        k = (tau + q0 + t + 2 * i) if case == 'MT' else (tau + q0 + t + 2 * i + 1)
        j = q0 + tau + i
        if not (30 <= k <= kmax and 3 <= j <= k - 2):
            continue
        sg = rng.choice(sigmas)
        C = tlib.Calc(k, j, tau, i, case, mode, sigma=sg)
        if not C.valid():
            continue
        names, norms, Q = tlib.reduced(C)
        if not names:
            continue
        n = len(names)
        conv = lambda x: x if isinstance(x, Fr) else (Fr(int(x.p), int(x.q)) if hasattr(x, 'p') else Fr(x))
        Qf = [[conv(Q[a][b]) for b in range(n)] for a in range(n)]
        kv = tlib.kernel_vectors(C, names)
        ann = all(sum(Qf[a][names.index(nm)] * Fr(cf) for nm, cf in v.items()) == 0 for _, v in kv for a in range(n))
        A = flint.fmpq_mat(n, n, [flint.fmpq(x.numerator, x.denominator) for row in Qf for x in row])
        rank = A.rank()
        null = len(kv)
        good = False
        if ann and rank == n - null:
            for drop in itertools.combinations(range(n), null):
                keep = [a for a in range(n) if a not in drop]
                sub = flint.fmpq_mat(len(keep), len(keep), [A[a, b] for a in keep for b in keep])
                if not keep or sub.det() != 0:
                    good = pd([[Qf[a][b] for b in keep] for a in keep])
                    break
        res = ann and rank == n - null and good
        ok &= res
        got += 1
        ndone += 1
        print(f'{case} {mode} {spec}: k={k} j={j} tau={tau} i={i} sigma={sg}: {n} vectors, {null} kernel vectors '
              f'annihilated {ann}, rank {rank}, PD complement {good}: {res} [{time.time()-T0:.0f}s]', flush=True)
print(f'{ndone} members of {len(fams)} certified families; {time.time()-T0:.0f}s')
print('PASS' if ok else 'FAIL')
