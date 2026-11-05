"""Negative controls (exact witnesses): with B replaced by (1 + 1/1000) B, every block of the certificate fails.

For a block with known kernel vectors kappa (rescaled coordinates), Y kappa = 0, so
kappa^T (A - (1 + e) B) kappa = -e kappa^T B kappa < 0 as soon as kappa^T B kappa > 0: checked exactly at a rational point.
For the blocks without kernel (j = 1, T = {0,1}; j = 2, T = {0,1,2}) the gap tends to 0 at the boundary of the orthant;
the script searches rational points with small weights and returns an explicit rational vector v with
v^T (A - (1 + e) B) v < 0, checked exactly.  PASS if every block has an exact witness."""
import time
import numpy as np
import flint
from k4lib import block_states, kernel_vectors, pt_builder, BLOCKS

t0 = time.time()
eps = flint.fmpq(1, 1000)
ok = True


def form(Mx, v):
    n = len(v)
    return sum(v[a] * Mx[a][b] * v[b] for a in range(n) for b in range(n) if v[a] != 0 and v[b] != 0)


def mats(Bd, S):
    A = [[Bd.A(a, yp, b, y) for (b, y) in S] for (a, yp) in S]
    B = [[Bd.B(a, yp, b, y) for (b, y) in S] for (a, yp) in S]
    return A, B


for (j, T) in BLOCKS:
    S = block_states(j, T)
    found = None
    if j == 1 and T == () or j == 1 and T == (0,) or j == 2 and T in ((), (0,), (0, 1)):
        s = [flint.fmpq(1), flint.fmpq(2, 3), flint.fmpq(5, 7), flint.fmpq(3, 11)]
        Bd = pt_builder(s)
        A, B = mats(Bd, S)
        kv = kernel_vectors(S, j, T, s)
        for v in kv:
            v = [flint.fmpq(x) for x in v]
            val = form([[A[a][b] - (1 + eps) * B[a][b] for b in range(len(S))] for a in range(len(S))], v)
            if val < 0:
                found = (s, v, val)
                break
    else:
        for e in range(1, 12):
            for shape in ((1, 1, 1), (1, 2, 3), (3, 1, 1), (1, 3, 1), (1, 1, 3)):
                for mode in ('small', 'large'):
                    x = flint.fmpq(1, 10 ** e)
                    if mode == 'small':
                        s = [flint.fmpq(1)] + [x * c for c in shape]
                    else:
                        s = [x] + [flint.fmpq(c) for c in shape]
                    Bd = pt_builder(s)
                    A, B = mats(Bd, S)
                    Ye = [[A[a][b] - (1 + eps) * B[a][b] for b in range(len(S))] for a in range(len(S))]
                    Yf = np.array([[float(z.p) / float(z.q) for z in row] for row in Ye])
                    ev, U = np.linalg.eigh(Yf)
                    if ev[0] < 0:
                        u = U[:, 0] / np.abs(U[:, 0]).max()
                        v = [flint.fmpq(int(round(c * 10 ** 6)), 10 ** 6) for c in u]
                        val = form(Ye, v)
                        if val < 0:
                            found = (s, v, val)
                            break
                if found:
                    break
            if found:
                break
    good = found is not None
    ok &= good
    if good:
        print(f'j = {j}, T = {T}: (1 + 1/1000) B fails at s = {[str(x) for x in found[0]]}: exact witness with '
              f'v^T Y_e v = {float(found[2].p) / float(found[2].q):.3e} < 0', flush=True)
    else:
        print(f'j = {j}, T = {T}: no witness found', flush=True)
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
