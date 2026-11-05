"""Check the closed form cf.TypeI against the exact Lemma B builder (lemmab.py) and against the direct float
construction (fastdirect.py)."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import numpy as np
from cf import TypeI
from lemmab import LemmaB
from fastdirect import Direct, typeI_states

ok = True
t0 = time.time()
cases = [((Fr(1), Fr(2), Fr(3)), [0]), ((Fr(1), Fr(2), Fr(3)), [1]), ((Fr(3, 2), Fr(1, 2), Fr(2), Fr(1)), [0]),
         ((Fr(3, 2), Fr(1, 2), Fr(2), Fr(1)), [0, 2]), ((Fr(3, 2), Fr(1, 2), Fr(2), Fr(1)), [1, 2, 3]),
         ((Fr(1), Fr(2), Fr(1, 3), Fr(5, 4), Fr(2, 3)), [0, 3]), ((Fr(1), Fr(2), Fr(1, 3), Fr(5, 4), Fr(2, 3)), [4]),
         ((Fr(1), Fr(2), Fr(1, 3), Fr(5, 4), Fr(2, 3)), [0, 1, 2])]
for w, T in cases:
    M = TypeI(w, T).build()
    Y, lab = LemmaB(w).typeI(T)
    err = max(abs(M[i][j] - 2 * Y[i][j]) for i in range(len(M)) for j in range(len(M)))
    good = err == 0
    ok &= good
    print(f"exact k={len(w)} T={T}: max |cf - 2 Y_lemmaB| = {err}  {'ok' if good else 'MISMATCH'}", flush=True)
rng = np.random.default_rng(5)
for k in (4, 5, 6):
    for trial in range(2):
        w = rng.random(k) + 0.2
        m = int(rng.integers(1, k - 1))
        T = sorted(rng.choice(k, m, replace=False).tolist())
        st, lab = typeI_states(k, T)
        Yd = Direct(w).matrix(st)
        M = np.array(TypeI(w, T, F=float).build())
        err = np.abs(M - 2 * Yd).max() / max(1, np.abs(M).max())
        good = err < 1e-10
        ok &= good
        print(f"float k={k} T={T}: rel err vs direct = {err:.2e}  {'ok' if good else 'MISMATCH'}", flush=True)
print(f"time {time.time() - t0:.1f}s")
print("PASS" if ok else "FAIL")
