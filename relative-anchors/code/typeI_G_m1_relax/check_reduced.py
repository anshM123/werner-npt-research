"""Exact check: the reduced system blockdiag D + U C U^T (reduced.pkl) equals the relaxed form Psi of psi_exact.py
at random rational weights, u = 3..8, after the change of coordinates (alpha, c, rho, zeta) -> (ah, ch, zh, rho):
ah = (Bm alpha)/w, ch = c/w, zh = zeta/(w_t w)."""
import random, time, sys
import flint
from fractions import Fraction as Fr
from reduced_c import ReducedC as Reduced, assemble, to_q
from psi_exact import PsiExact

Q = flint.fmpq
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 2024)
t0 = time.time()
R = Reduced()
print(f'loaded reduced system: J = {len(R.gens)} global directions (rho and {len(R.gens) - 1} functionals), '
      f'{len(R.Cterms)} C-terms, {len(R.loc)} D-terms ({time.time() - t0:.1f} s)', flush=True)
ok = True
for u in range(3, 9):
    for trial in range(2):
        wt = Fr(random.randint(1, 40), random.randint(1, 40))
        wU = [Fr(random.randint(1, 40), random.randint(1, 40)) for _ in range(u)]
        num = R.numeric(wt, wU)
        Psi_r, Dm, Um, Cm = assemble(num)
        P = PsiExact(Q(wt.numerator, wt.denominator), [Q(x.numerator, x.denominator) for x in wU])
        Pe, iy = P.matrix()
        n_pe = Pe.nrows()
        n_r = 3 * u + 1
        T = flint.fmpq_mat(n_r, n_pe)
        w = [Q(x.numerator, x.denominator) for x in wU]
        wtq = Q(wt.numerator, wt.denominator)
        for i in range(u):
            for k in range(u - 1):
                T[i, k] = P.Bm[i, k] / w[i]
            T[u + i, (u - 1) + i] = 1 / w[i]
            T[2 * u + i, (u - 1) + u + 1 + i] = 1 / (wtq * w[i])
        T[3 * u, (u - 1) + u] = Q(1)
        diff = T.transpose() * Psi_r * T - Pe
        good = all(diff[i, j] == 0 for i in range(n_pe) for j in range(n_pe))
        ok &= good
        print(f'u = {u}, trial {trial}: T^T (D + U C U^T) T == Psi exactly: {good}', flush=True)
print(f'total {time.time() - t0:.1f} s')
print('PASS' if ok else 'FAIL')
