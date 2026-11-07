"""Exact verification (rational arithmetic, coupled basis) of a word-span certificate at a = 3/4, tau' = 1.
usage: python verify_exact.py n greedy
       python verify_exact.py n FILE.pkl [maxden]     (pickle with 'free' pairs and float 'y'; rationalised)
Checks: (1) formal identity Q + Gamma(Q) = K' - I as polynomials in p, q, r (exact);
        (2) for every spin L the rescaled block Q'' is positive semidefinite (exact LDL^T); total kernel = n + 1."""
import sys, time, math, pickle
from fractions import Fraction
from exact_cpl import ExactL, q_block, ldl_psd
from wsdp import poly_W, multinom, pairs, add, K_poly, ident_poly

def certificate_poly(n, yfree):
    """Q = sum_{s>t} C[y W(s,t) + (1-y) W(t,s)] + 1/2 sum_{s=t>=1} C W(s,s);  y = 1 unless given in yfree."""
    P = {}
    for (s, t) in pairs(n):
        y = Fraction(yfree.get((s, t), 1))
        C = multinom(n, n - s - t, s, t)
        P = add(P, poly_W(n, s, t), C * y)
        P = add(P, poly_W(n, t, s), C * (1 - y))
    for s in range(1, n // 2 + 1):
        P = add(P, poly_W(n, s, s), Fraction(multinom(n, n - 2 * s, s, s), 2))
    return {k: Fraction(v) for k, v in P.items() if v != 0}

def main():
    n = int(sys.argv[1])
    if sys.argv[2] == 'greedy':
        yfree = {}
    else:
        d = pickle.load(open(sys.argv[2], 'rb'))
        maxden = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
        yfree = {p: Fraction(float(v)).limit_denominator(maxden) for p, v in zip(d['free'], d['y'])}
    t0 = time.time()
    P = certificate_poly(n, yfree)
    G = {(a, c, b): v for (a, b, c), v in P.items()}
    lhs = add(P, G)
    rhs = add(K_poly(n), ident_poly(n), -1)
    assert all(v == 0 for v in add(lhs, rhs, -1).values()), 'identity fails'
    print(f"n={n}: certificate {'greedy' if not yfree else 'non-greedy pairs ' + str({p: str(v) for p, v in yfree.items() if v != 1})}", flush=True)
    print(f"   identity Q + Gamma(Q) = K' - I holds exactly", flush=True)
    allok = True; ker = 0
    for L in range(0, 2 * n + 1):
        EL = ExactL(n, L)
        Q = q_block(EL, P)
        ok, rank, zeros = ldl_psd(Q)
        allok &= ok; ker += zeros
        if not ok:
            print(f"   L={L}: NOT positive semidefinite (size {Q.nrows()})", flush=True)
    print(f"   RESULT n={n}: Q positive semidefinite (exact) = {allok}; kernel dimension (sum over L of zero pivots, "
          f"meaningful only if PSD) = {ker} (n + 1 = {n + 1}); {time.time() - t0:.1f}s", flush=True)

if __name__ == '__main__':
    main()
