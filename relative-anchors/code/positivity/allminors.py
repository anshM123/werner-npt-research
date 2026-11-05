"""Test 3 (exact): every principal minor of a gauge-fixed block, with the denominator factors pulled out.

Nt = L * den * M_g (integer coefficients, den = lcm of the entry denominators, L a positive integer).  For a subset S
in increasing order, Sylvester's identity gives the bordered determinants
    B_P[a][b] = det Nt[P + a, P + b],   B_{P+i}[a][b] = (B_P[i][i] B_P[a][b] - B_P[a][i] B_P[i][b]) / det Nt[P, P],
with exact division, and det Nt[P + i, P + i] = B_P[i][i].  A depth-first search over the subsets (children of P are
P + i, i > max P) computes every principal minor once.  Entries are stored as core * prod f_i^{x_i} with the
irreducible factors f_i of den divided out exactly (as in ../middle_k4/bareiss_cert.py).  Every minor is
    det M_g[S, S] = core_S * prod_i f_i^{x_iS - |S| E_i} / (cst L)^{|S|},
so it is subtraction-free if cst > 0, every f_i has positive coefficients and core_S has positive coefficients.
Each minor is also checked against a direct exact determinant of Nt[S, S] at a random integer point.
Usage: python allminors.py k [blockname|all] [maxsize] [random:N:seed | sub:N:m:seed]
(random: N leading-minor chains in random orders; sub: N random m-subsets, each in a random order, with its
m nested leading minors)"""
import sys, time, random, itertools
from math import lcm
from fractions import Fraction as Fr
import flint
from psdtools import load_blocks, poly_matrix



def detF(A):
    A = [row[:] for row in A]
    n = len(A)
    d = Fr(1)
    for p in range(n):
        piv = next((i for i in range(p, n) if A[i][p] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != p:
            A[p], A[piv] = A[piv], A[p]
            d = -d
        d *= A[p][p]
        for i in range(p + 1, n):
            f = A[i][p] / A[p][p]
            for l in range(p, n):
                A[i][l] -= f * A[p][l]
    return d


def run_block(blk, randspec=None, quiet=False):
    n = len(blk['Mg'])
    names = blk['names']
    ctxZ = flint.fmpz_mpoly_ctx.get(tuple(names), 'lex')
    den, Nt = poly_matrix(blk)
    L = 1
    for row in Nt:
        for x in row:
            for c in x.coeffs():
                L = lcm(L, int(c.q))
    NtZ = [[ctxZ.from_dict({m: int(c * L) for m, c in x.to_dict().items()}) for x in row] for row in Nt]
    denZ = ctxZ.from_dict({m: int(c * lcm(*[int(v.q) for v in den.coeffs()])) for m, c in den.to_dict().items()})
    cst, fac = denZ.factor()
    FACS = [f for f, e in fac]
    nf = len(FACS)
    facpos = cst > 0 and all(all(c > 0 for c in f.coeffs()) for f in FACS)
    ZERO = ctxZ.from_dict({})
    ONE = ctxZ.from_dict({(0,) * len(names): 1})

    def strip(p):
        ex = [0] * nf
        if p == 0:
            return p, ex
        for i, f in enumerate(FACS):
            while True:
                try:
                    q = p / f
                except Exception:
                    break
                p = q
                ex[i] += 1
        return p, ex

    rng = random.Random(7)
    pt = [rng.randint(2, 30) for _ in names]
    Nv = [[Fr(int(NtZ[a][c](*pt))) for c in range(n)] for a in range(n)]
    fval = [Fr(int(f(*pt))) for f in FACS]
    stats = dict(n=0, pos=0, mixed=0, badcheck=0, maxterms=0, maxdeg=0)
    mixed_examples = []

    def record(S, core, ex):
        stats['n'] += 1
        val = Fr(int(core(*pt)))
        for fv, e in zip(fval, ex):
            val *= fv ** e
        if val != detF([[Nv[a][c] for c in S] for a in S]):
            stats['badcheck'] += 1
        cs = core.coeffs()
        if all(c > 0 for c in cs):
            stats['pos'] += 1
        else:
            stats['mixed'] += 1
            if len(mixed_examples) < 5:
                mixed_examples.append((tuple(S), core.total_degree(), len(core), sum(1 for c in cs if c < 0)))
        stats['maxterms'] = max(stats['maxterms'], len(core))
        stats['maxdeg'] = max(stats['maxdeg'], core.total_degree())

    def step(Bp, prev, i, idx):
        """Bp: dict (a, b) -> (core, ex) for a <= b in idx (indices > max P); returns B_{P+i} on idx' = [x > i]."""
        cp, ep = Bp[(i, i)]
        cprev, eprev = prev
        idx2 = [x for x in idx if x > i]
        out = {}
        for ai, a in enumerate(idx2):
            ca_i, ea_i = Bp[(min(a, i), max(a, i))]
            for b in idx2[ai:]:
                cb_i, eb_i = Bp[(min(b, i), max(b, i))]
                cab, eab = Bp[(a, b)]
                e1 = [x + y for x, y in zip(ep, eab)]
                e2 = [x + y for x, y in zip(ea_i, eb_i)]
                if cab == 0:
                    g = e2
                elif ca_i == 0 or cb_i == 0:
                    g = e1
                else:
                    g = [min(x, y) for x, y in zip(e1, e2)]
                t = ZERO
                if cab != 0:
                    m1 = cp * cab
                    for f, x, gg in zip(FACS, e1, g):
                        if x > gg:
                            m1 = m1 * f ** (x - gg)
                    t = t + m1
                if ca_i != 0 and cb_i != 0:
                    m2 = ca_i * cb_i
                    for f, x, gg in zip(FACS, e2, g):
                        if x > gg:
                            m2 = m2 * f ** (x - gg)
                    t = t - m2
                if t == 0:
                    out[(a, b)] = (ZERO, [0] * nf)
                    continue
                q = t / cprev
                ex = [gg - x for gg, x in zip(g, eprev)]
                for kk, f in enumerate(FACS):
                    while ex[kk] < 0:
                        q = q / f
                        ex[kk] += 1
                c2, e3 = strip(q)
                out[(a, b)] = (c2, [x + y for x, y in zip(ex, e3)])
        return out, idx2

    B0 = {}
    for a in range(n):
        for b in range(a, n):
            B0[(a, b)] = strip(NtZ[a][b])

    def dfs(P, Bp, prev, idx):
        for i in idx:
            S = P + [i]
            core, ex = Bp[(i, i)]
            record(S, core, ex)
            if len(S) < maxdepth and any(x > i for x in idx):
                B2, idx2 = step(Bp, prev, i, idx)
                dfs(S, B2, (core, ex), idx2)

    maxdepth = n
    t0 = time.time()
    if randspec is None:
        dfs([], B0, (ONE, [0] * nf), list(range(n)))
        mode = 'all subsets'
    else:
        parts = randspec.split(':')
        if parts[0] == 'sub':
            _, N, msz, seed = parts
            msz = min(int(msz), n)
        else:
            _, N, seed = parts
            msz = n
        rr = random.Random(int(seed))
        mode = f'{N} random chains of length {msz} (seed {seed})'
        for _ in range(int(N)):
            perm = list(range(n))
            rr.shuffle(perm)
            perm = perm[:msz]
            n_save = n
            # a random chain: leading minors in a random order = nested subsets; computed by relabelling
            Bp = {}
            for a in range(msz):
                for b in range(a, msz):
                    pa, pb = perm[a], perm[b]
                    Bp[(a, b)] = strip(NtZ[pa][pb])
            prev = (ONE, [0] * nf)
            idx = list(range(msz))
            for r in range(msz):
                core, ex = Bp[(r, r)]
                record([perm[x] for x in range(r + 1)], core, ex)
                if r < msz - 1:
                    Bp2, idx2 = step(Bp, prev, r, idx)
                    prev = (core, ex)
                    Bp, idx = Bp2, idx2
    if quiet:
        return stats, facpos
    print(f'{blk["name"]} ({blk["kind"]}), size {n}, {mode}: {stats["n"]} principal minors, cores with positive '
          f'coefficients {stats["pos"]}, mixed {stats["mixed"]}, point-check failures {stats["badcheck"]}; den factors '
          f'positive {facpos}; largest core degree {stats["maxdeg"]}, {stats["maxterms"]} terms  [{time.time() - t0:.1f}s]',
          flush=True)
    for ex in mixed_examples:
        print(f'    mixed core: subset {[blk["labels"][x] for x in ex[0]]}, degree {ex[1]}, {ex[2]} terms, {ex[3]} negative')
    return stats, facpos


if __name__ == '__main__':
  k = int(sys.argv[1])
  which = sys.argv[2] if len(sys.argv) > 2 else 'all'
  maxsize = int(sys.argv[3]) if len(sys.argv) > 3 else 99
  randspec = sys.argv[4] if len(sys.argv) > 4 else None
  t00 = time.time()
  ok = True
  tot = dict(n=0, pos=0, mixed=0, bad=0)
  for blk in load_blocks(k):
    if which != 'all' and blk['name'] != which:
        continue
    if len(blk['Mg']) > maxsize:
        continue
    st, fp = run_block(blk, randspec)
    tot['n'] += st['n']
    tot['pos'] += st['pos']
    tot['mixed'] += st['mixed']
    tot['bad'] += st['badcheck']
    ok &= st['badcheck'] == 0 and fp
  print(f'TOTAL: {tot["n"]} principal minors, positive cores {tot["pos"]}, mixed {tot["mixed"]}, point-check failures '
        f'{tot["bad"]}  [{time.time() - t00:.1f}s]')
  print('PASS' if ok else 'FAIL')
