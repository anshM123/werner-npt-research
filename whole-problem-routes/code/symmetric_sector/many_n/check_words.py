"""Brute-force check of words.py against the earlier monomial-basis code (bos.py, ybasis.py/zbasis.py, o3red.py):
for n = 1..5, every word M(k,s,t) (one arrangement of I^k F^s Omega^t, compressed) has the same spectrum on each spin-L
block, Gamma(M(k,s,t)) = M(k,t,s) exactly in the monomial basis, and the blocks of K_n = Pi C_{3/4}^{(x)n} Pi agree."""
import sys, math
import numpy as np
from zbasis import one_copy_ops, transfer
from bos import build_K
from o3red import O3Red
from words import level_data, word_block

def build_M(n):
    ops = one_copy_ops(3)
    lst = [ops['I'], ops['F'], ops['Om']]
    Y = {(0, 0, 0): np.array([[1.0]])}
    for m in range(1, n + 1):
        Ynew = {}
        for e1 in range(m + 1):
            for e2 in range(m + 1 - e1):
                e = (e1, e2, m - e1 - e2)
                for tp in range(3):
                    if e[tp] > 0:
                        pr = list(e); pr[tp] -= 1
                        Ynew[e] = transfer(Y[tuple(pr)], m, lst[tp], 3)
                        break
        Y = Ynew
    return Y

worst = 0.0
for n in range(1, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 5):
    R = O3Red(n)
    Mbf = build_M(n)
    N = R.N
    G = lambda X: X.reshape(N, N, N, N).transpose(0, 3, 2, 1).reshape(N * N, N * N)
    for (k, s, t), X in Mbf.items():
        worst = max(worst, np.abs(G(X) - Mbf[(k, t, s)]).max())
    Kbf = build_K(n, 0.75) * 4 ** n
    for L in R.Lvals:
        ld = level_data(n, L)
        Kc = 0
        for (k, s, t), X in Mbf.items():
            Wc = word_block(ld, n, k, s, t)
            ev1 = np.linalg.eigvalsh(Wc); ev2 = np.linalg.eigvalsh(R.block(X, L))
            worst = max(worst, np.abs(ev1 - ev2).max())
            Kc = Kc + math.factorial(n) // (math.factorial(k) * math.factorial(s) * math.factorial(t)) * 3 ** k * (-1) ** (s + t) * Wc
        ev1 = np.linalg.eigvalsh(Kc); ev2 = np.linalg.eigvalsh(R.block(Kbf, L))
        worst = max(worst, np.abs(ev1 - ev2).max() / max(1, np.abs(ev2).max()))
    print(f"n={n}: words, Gamma and K_n blocks agree with brute force; worst deviation so far {worst:.2e}", flush=True)
print("PASS" if worst < 1e-8 else "FAIL")
