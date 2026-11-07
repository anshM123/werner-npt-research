"""Checks of cpl.py: spectra of E_ab E_ba and (a.b)^T(a.b) on every spin-L block against the closed forms
(m-j)(m-j+1) [gl(2) Casimir, GL(3) type (2m-j, j)] and tau(2n - tau + 2) [su(1,1) Casimir]; dimension counts."""
import sys
import numpy as np
from collections import Counter
from cpl import basis, ab, Eba, mult

def dimGL3(l1, l2):
    return (l1 - l2 + 1) * (l1 + 2) * (l2 + 1) // 2

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
worst = 0.0
for m in range(1, nmax + 1):
    # gl(2) Casimir check
    cnt = Counter()
    for L in range(0, 2 * m + 1):
        E = Eba(m, L)
        C = E.T @ E
        ev = np.linalg.eigvalsh(C)
        allowed = np.array([(m - j) * (m - j + 1) for j in range(m + 1)], float)
        for e in ev:
            k = int(np.argmin(np.abs(allowed - e)))
            worst = max(worst, abs(allowed[k] - e))
            cnt[k] += 2 * L + 1
    for j in range(m + 1):
        assert cnt[j] == dimGL3(2 * m - j, j), (m, j, cnt[j], dimGL3(2 * m - j, j))
    # su(1,1) Casimir check
    for L in range(0, 2 * m + 1):
        A = ab(m, L)
        ev = np.sort(np.linalg.eigvalsh(A.T @ A))
        exp = []
        for tau in range(0, m + 1):
            mult_tau = len(basis(m - tau, m - tau, L)) - len(basis(m - tau - 1, m - tau - 1, L))
            exp += [tau * (2 * m - tau + 2)] * mult_tau
        exp = np.sort(np.array(exp, float))
        assert len(exp) == len(ev)
        worst = max(worst, np.abs(exp - ev).max() if len(ev) else 0)
    print(f"m={m}: gl(2) and su(1,1) Casimir spectra OK; block sizes {[mult(m, L) for L in range(2*m+1)]}; worst dev {worst:.2e}", flush=True)
print("PASS" if worst < 1e-8 else "FAIL")
