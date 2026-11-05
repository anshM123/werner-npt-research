"""Coverage: every block of each family and every weight vector with at most two distinct weights reduce to a
certified configuration.

A weight vector with at most two distinct squared weights is given by a set A of copies (squared weight alpha) and
its complement (squared weight beta), alpha, beta > 0 (alpha = beta, A = all copies or A empty: one level).  For a
block labelled by T (|T| = m) put a = |T n A|, b = |U n A| (U the other copies), and let pi be the permutation that
maps T n A, T - A, U n A, U - A in increasing order onto 0..a-1, a..m-1, m..m+b-1, m+b..k-1.  Then pi maps the states
of the block labelled by T onto those labelled by T0 = {0..m-1}, and, the rescaled block being homogeneous of degree
d = j + 1 - m in the squared weights,
    Yr(T, s)[X1, X2] = alpha^d Yr(T0, configuration (a, b) at q = beta/alpha)[pi X1, pi X2].
So the block is >= 0 with exactly the known kernel if configuration (a, b) is certified for every q > 0.  If only
the exchanged configuration (m - a, k - m - b) is certified, the same holds with A replaced by its complement and
alpha, beta exchanged (ratio alpha/beta, factor beta^d).
This script checks, for every family, every T and every A:
 - that the configuration (a, b), or its exchange, was certified (logs/run_cert_k<k>_<family>_<a>_<b>.log ends with
   PASS);
 - the identity above exactly at alpha = 4/9, beta = 25/16, with the left side built by j1lib.py (Lemma B) and the
   right side from the saved certified block (logs/cert_k<k>_<family>_<a>_<b>.pkl): for every pattern, or with
   sample=N for every N-th pattern.
python coverage.py FAMILIES [k=5] [sample=N]
Prints PASS or FAIL on the last line."""
import sys, time, itertools, os
from fractions import Fraction as Fr
import k5lib as K
from saved import load, block_at, cert_path

t0 = time.time()
args = [a for a in sys.argv[1:] if '=' not in a]
opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
k = int(opts.get('k', 5))
fams = args[0].split(',')
ok = True
alpha, beta = Fr(4, 9), Fr(25, 16)
wa, wb = K.rsqrt(alpha), K.rsqrt(beta)


def perm_string(y, pi):
    z = [0] * len(y)
    for t, x in enumerate(y):
        z[pi[t]] = x
    return tuple(z)


for f in fams:
    j, m = K.fam(f)
    d = j + 1 - m
    certified = set()
    for (a, b) in K.configs(k, m):
        lg = f'logs/run_cert_k{k}_{f}_{a}_{b}.log'
        if os.path.exists(lg):
            lines = [x.strip() for x in open(lg).read().splitlines() if x.strip()]
            if lines and lines[-1] == 'PASS':
                certified.add((a, b))
    cache = {}
    npat = nid = nexch = 0
    used = set()
    sample = opts.get('sample', 'all')
    for T in itertools.combinations(range(k), m):
        U = [t for t in range(k) if t not in T]
        for r in range(k + 1):
            for A in itertools.combinations(range(k), r):
                A = set(A)
                la, lb = alpha, beta            # squared weight on A, on the complement
                a, b = len([t for t in T if t in A]), len([t for t in U if t in A])
                if (a, b) not in certified:
                    # level exchange: the complement of A is the level of squared weight 1, at ratio alpha/beta
                    A = set(range(k)) - A
                    la, lb = beta, alpha
                    a, b = len([t for t in T if t in A]), len([t for t in U if t in A])
                    nexch += 1
                if (a, b) not in certified:
                    ok = False
                    print(f'  {f} T = {T} A = {sorted(A)}: configuration {(a, b)} NOT certified (nor its exchange)')
                    continue
                TA = [t for t in T if t in A]
                TB = [t for t in T if t not in A]
                UA = [t for t in U if t in A]
                UB = [t for t in U if t not in A]
                pi = {}
                for i, t in enumerate(TA + TB + UA + UB):
                    pi[t] = i
                assert sorted(pi.values()) == list(range(k))
                used.add((a, b))
                npat += 1
                # exact identity, for every pattern (sample = all) or for every sample-th pattern
                if sample != 'all' and (npat - 1) % int(sample) != 0:
                    continue
                w = [K.rsqrt(la) if t in A else K.rsqrt(lb) for t in range(k)]
                st, M = K.build_j1(k, j, T, w)
                key = (a, b, lb / la)
                if key not in cache:
                    data = load(cert_path(f, a, b, k))
                    cache[key] = (data['states'], block_at(data, lb / la))
                st0, M0 = cache[key]
                idx0 = {x: i for i, x in enumerate(st0)}
                img = [idx0.get((perm_string(a_, pi), perm_string(y, pi))) for (a_, y) in st]
                good = None not in img and sorted(img) == list(range(len(st0)))
                if good:
                    good = all(M[i][l] == la ** d * M0[img[i]][img[l]] for i in range(len(st)) for l in range(len(st)))
                ok &= good
                nid += 1
                if not good:
                    print(f'  {f} T = {T} A = {sorted(A)}: identity FAILS')
    print(f'{f}: {npat} (T, A) patterns ({nexch} through the level exchange), each reduces to a certified '
          f'configuration: {ok}; exact identity checked for {nid} of them; configurations used {sorted(used)} '
          f'[{time.time() - t0:.1f}s]', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
