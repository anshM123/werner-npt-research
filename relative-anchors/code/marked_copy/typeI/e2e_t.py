"""End to end at points, exact, without the reduction: whole type I blocks with a letter 1 from the Lemma B builder.

python e2e_t.py K W0S [JS]
At k = K, for every j in JS (default 3..K-1), every tau = 1..j and both cases (the marked copy 0 in T = {0..tau-1}, or
in U with T = {1..tau}), and every w0 in W0S (squared weight sigma = w0^2 at copy 0, 1 elsewhere):
  - the rescaled block of explicit_t.block (fastj1.FastJ1, Lemma B, exact over Q);
  - the known kernel vectors p (x) y (one per column y of the letter-1 rows) are annihilated exactly;
  - the block without a set of pivot coordinates of them (exact RREF) is positive definite by the exact test
    crosscheck.pd_cert (integer congruence, strict diagonal dominance in exact integer arithmetic);
so the block is >= 0 with kernel exactly span(p (x) y) at these points.  Negative control (B -> (1 + 1/1000) B, same
blocks, first w0): <kappa, Y~ kappa> < 0 exactly for the first known kernel vector kappa, so the scaled block is not
positive semidefinite."""
import sys, time
from fractions import Fraction as Fr
from explicit_t import block, block_states, known_kernel, layout
from crosscheck import pd_cert, pivots_of


def main():
    k = int(sys.argv[1])
    w0s = [Fr(x) for x in sys.argv[2].split(',')]
    js = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else list(range(3, k))
    T0 = time.time()
    ok = True
    nb = 0
    for j in js:
        for tau in range(1, j + 1):
            for case in ('MT', 'MU'):
                if case == 'MU' and tau > k - 1:
                    continue
                T = layout(k, tau, case, 0)[0]
                for wn, w0 in enumerate(w0s):
                    st, Y = block(k, j, T, w0)
                    n = len(st)
                    kv = known_kernel(k, j, T, st, w0)
                    kern = all(sum(Y[x][c] * v for c, v in vec.items()) == 0 for vec in kv for x in range(n))
                    vecs = [[vec.get(c, Fr(0)) for c in range(n)] for vec in kv]
                    piv = set(pivots_of(vecs)) if vecs else set()
                    keep = [x for x in range(n) if x not in piv]
                    proved, info = pd_cert([[Y[a][b] for b in keep] for a in keep])
                    nb += 1
                    good = kern and proved
                    neg = None
                    if wn == 0 and kv:
                        st2, Ys = block(k, j, T, w0, bscale=Fr(1001, 1000))
                        vec = kv[0]
                        val = sum(c1 * Ys[a][b] * c2 for a, c1 in vec.items() for b, c2 in vec.items())
                        neg = val < 0
                        good = good and neg
                    ok &= good
                    print(f'k={k} j={j} tau={tau} {case} sigma={w0 * w0}: {n} states, {len(kv)} kernel vectors '
                          f'annihilated {kern}, rest PD {proved} ({info}); scaled control negative {neg} '
                          f'[{time.time() - T0:.1f}s]', flush=True)
    print(f'{nb} blocks')
    print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')


if __name__ == '__main__':
    main()
