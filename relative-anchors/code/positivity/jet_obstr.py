"""Test 7 (exact negative statement): no fixed-vector representation with subtraction-free weights.

Suppose M_g = sum_e c_e x_e x_e^T with fixed real vectors x_e and weights c_e that are ratios of polynomials with
positive coefficients.  With Q = prod_e (denominator of c_e), Q * M_g has every monomial coefficient matrix PSD, and
Q = q * den with q a polynomial (den = lcm of the entry denominators) whose top coefficient in every generic direction
is positive (Q and den have positive coefficients).  Let C_m be the coefficient matrices of den * M_g and l a direction
for which the maximum of l.m over the support is attained at a single monomial v and the second largest value at a
single monomial m1.  The coefficient matrix of Q * M_g at u0 + m1 (u0 the top monomial of q) is
q_{u0} C_{m1} + (coefficient of q) * C_v, so
  (first order)   C_v must be PSD, and
  (second order)  C_{m1} restricted to ker C_v must be PSD.
A violation, with the direction l and an exact witness, excludes every such representation.
Usage: python jet_obstr.py k [maxsize]"""
import sys, itertools, time
import flint
from fractions import Fraction as Fr
from psdtools import load_blocks, poly_matrix, coeff_mats, psd_exact

k = int(sys.argv[1])
maxsize = int(sys.argv[2]) if len(sys.argv) > 2 else 99
t0 = time.time()
found_blocks = []
tested = 0
for blk in load_blocks(k):
    n = len(blk['Mg'])
    if n > maxsize:
        continue
    den, Nt = poly_matrix(blk)
    C = coeff_mats(Nt)
    if all(psd_exact(C[m])[0] for m in C):
        continue
    tested += 1
    supp = list(C)
    nv = len(supp[0])
    rng = 4 if nv <= 3 else 3
    cert = None
    for l in itertools.product(range(-rng, rng + 1), repeat=nv):
        if not any(l):
            continue
        vals = sorted(((sum(x * y for x, y in zip(l, m)), m) for m in supp), reverse=True)
        if len(vals) < 2 or vals[0][0] == vals[1][0] or (len(vals) > 2 and vals[1][0] == vals[2][0]):
            continue
        v, m1 = vals[0][1], vals[1][1]
        okv, wv = psd_exact(C[v])
        if not okv:
            cert = ('first order', l, v, None)
            break
        # kernel of C_v
        Cv = flint.fmpq_mat(n, n, [flint.fmpq(x.numerator, x.denominator) for row in C[v] for x in row])
        Rr, rank = Cv.rref()
        piv = []
        r_ = 0
        for c in range(n):
            if r_ < rank and Rr[r_, c] != 0:
                piv.append(c)
                r_ += 1
        free = [c for c in range(n) if c not in piv]
        if not free:
            continue
        K = []
        for f in free:
            vec = [Fr(0)] * n
            vec[f] = Fr(1)
            for i, pc in enumerate(piv):
                x = -Rr[i, f]
                vec[pc] = Fr(int(x.p), int(x.q))
            K.append(vec)
        Cm = C[m1]
        R = [[sum(K[i][a] * Cm[a][b] * K[j][b] for a in range(n) for b in range(n)) for j in range(len(K))]
             for i in range(len(K))]
        okr, wr = psd_exact(R)
        if not okr:
            x = [sum(wr[i] * K[i][a] for i in range(len(K))) for a in range(n)]
            val = sum(x[a] * Cm[a][b] * x[b] for a in range(n) for b in range(n))
            assert val < 0 and all(sum(C[v][a][b] * x[b] for b in range(n)) == 0 for a in range(n))
            cert = ('second order', l, v, m1, [(blk['labels'][a], str(x[a])) for a in range(n) if x[a] != 0])
            break
    if cert:
        found_blocks.append(blk['name'])
        print(f'{blk["name"]} ({blk["kind"]}), size {n}: {cert[0]} obstruction in direction {cert[1]}: top monomial '
              f'{cert[2]}, next {cert[3] if len(cert) > 3 else ""}; witness in ker C_v with x^T C_m1 x < 0: '
              f'{cert[4] if len(cert) > 4 else ""}  [{time.time() - t0:.1f}s]', flush=True)
    else:
        print(f'{blk["name"]} ({blk["kind"]}), size {n}: no first- or second-order obstruction over the directions '
              f'tested  [{time.time() - t0:.1f}s]', flush=True)
print(f'RESULT k = {k}: {len(found_blocks)} of {tested} blocks without a q = 1 representation have an exact jet '
      f'obstruction: {found_blocks}')
print('PASS')
