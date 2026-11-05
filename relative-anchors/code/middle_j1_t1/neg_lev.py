"""Negative controls (B -> (1 + 1/1000) B) for the level certificates.

The B-part of every entry (every term with an E^{-1}) is multiplied by 1001/1000.  For equal U-weights (l = 1), two
U-weight levels (l = 2) and three levels (two single U-copies and two equal ones, l = 3) we exhibit explicit rational
weight vectors at which the scaled gauge-fixed block (state P_t deleted) is not positive definite, while the true one
is:
 (1) from the symbolic forms Ginv, Gstd of lev.py, evaluated exactly at the point: some leading minor of the scaled
     invariant or standard form is negative, and all leading minors of the true forms are positive;
 (2) independently, from j1lib (Lemma B, exact): the whole gauge-fixed block with bscale 1001/1000 is not positive
     definite (exact LDL), while with bscale 1 it is.
(The certificate itself fails for the scaled B-part: log of cert_lev.py G 1 2 1001/1000.)
Usage: python neg_lev.py"""
import time
from fractions import Fraction as Fr
import flint
from lev import build
from j1lib import J1, block_states

t0 = time.time()


def ev(f, vals):
    def pr(m):
        out = Fr(1)
        for v, e in zip(vals, m):
            out *= v ** int(e)
        return out
    num = sum((Fr(int(flint.fmpq(c).p), int(flint.fmpq(c).q)) * pr(m) for m, c in f.n.to_dict().items()), Fr(0))
    den = sum((Fr(int(flint.fmpq(c).p), int(flint.fmpq(c).q)) * pr(m) for m, c in f.d.to_dict().items()), Fr(0))
    return num / den


def ldl_pd(A):
    A = [row[:] for row in A]
    n = len(A)
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                for j in range(c + 1, n):
                    A[i][j] -= f * A[c][j]
    return True


def lead_minors(A):
    """leading principal minors over Q (Gaussian elimination without pivoting; stops at a zero pivot)."""
    A = [row[:] for row in A]
    n = len(A)
    out = []
    d = Fr(1)
    for c in range(n):
        if A[c][c] == 0:
            out.append(Fr(0))
            break
        d *= A[c][c]
        out.append(d)
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            for j in range(c, n):
                A[i][j] -= f * A[c][j]
    return out


found = []
for l, offs, us, pts in ((1, [2], [2], [(Fr(1), [Fr(1, 10)]), (Fr(1), [Fr(1, 30)]), (Fr(1), [Fr(1, 100)])]),
                         (1, [2], [5], [(Fr(1), [Fr(1, 30)]), (Fr(1), [Fr(1, 100)])]),
                         (2, [1, 1], [2, 2], [(Fr(1), [Fr(1, 30), Fr(1, 60)]), (Fr(1), [Fr(1, 100), Fr(1, 300)])]),
                         (2, [1, 1], [1, 3], [(Fr(1), [Fr(1, 50), Fr(1, 100)]), (Fr(1), [Fr(1, 100), Fr(1, 200)])]),
                         (3, [1, 1, 1], [1, 1, 2], [(Fr(1), [Fr(1, 50), Fr(1, 70), Fr(1, 100)]),
                                                    (Fr(1), [Fr(1, 100), Fr(1, 150), Fr(1, 200)])])):
    S1 = build(l, offs, 1, want_phi=False)
    Sb = build(l, offs, Fr(1001, 1000), want_phi=False)
    for wt, wl in pts:
        vals = [wt * wt] + [v * v for v in wl] + [Fr(us[i] - offs[i]) for i in range(l)]
        forms1 = [S1['Ginv']] + [S1['Gstd'][i] for i in range(l) if us[i] >= 2]
        formsb = [Sb['Ginv']] + [Sb['Gstd'][i] for i in range(l) if us[i] >= 2]
        m1 = [lead_minors([[ev(f, vals) for f in row] for row in A]) for A in forms1]
        mb = [lead_minors([[ev(f, vals) for f in row] for row in A]) for A in formsb]
        true_pos = all(v > 0 for ms in m1 for v in ms)
        scaled_neg = any(v <= 0 for ms in mb for v in ms)
        w = [wt] + [wl[i] for i in range(l) for _ in range(us[i])]
        k = len(w)
        st = block_states(k, (0,))
        keep = [i for i, (a, y) in enumerate(st) if not (a is not None and a[0] == 1 and a[1] == 0)]
        Mt = J1(w).matrix(st)
        Ms = J1(w, bscale=Fr(1001, 1000)).matrix(st)
        pd_true = ldl_pd([[Mt[i][j] for j in keep] for i in keep])
        pd_scaled = ldl_pd([[Ms[i][j] for j in keep] for i in keep])
        good = true_pos and scaled_neg and pd_true and not pd_scaled
        print(f'l = {l}, u = {us}, w_t = {wt}, U-weights {[str(v) for v in wl]}: symbolic forms: true minors positive '
              f'{true_pos}, scaled minors with a nonpositive one {scaled_neg}; j1lib gauge-fixed block PD: true {pd_true}, '
              f'scaled {pd_scaled} -> witness {good} [{time.time() - t0:.1f}s]', flush=True)
        if good:
            found.append((l, tuple(us)))
ok = all(any(f[0] == lv for f in found) for lv in (1, 2, 3))
print(f'witnesses: {found}')
print('PASS' if ok else 'FAIL')
