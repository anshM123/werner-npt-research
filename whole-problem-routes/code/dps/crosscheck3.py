"""Floating-point end-to-end check of the conventions in verify_independent.py (three copies).
The left side is computed from the definition of the map P_a only (no diagrams):  for random complex x, y in C^27,
    LHS = <y| P_a^{(x)3}(|conj x><conj x|) |y> |y|^2 - t |x|^2 |y|^4,
and the right side from the exact certificate operators of verify_independent.py:
    RHS = <psi|Z|psi> + <conj(x) y y|QA'|conj(x) y y> + <x conj(y) y|QB|x conj(y) y>,   psi = x (x) y (x) y.
They must agree to rounding, and each term of RHS must be >= 0.  Also checks <x(x)y|J|x(x)y> = LHS-part with J the
Choi matrix of P_a^{(x)3} built entry by entry.
usage: python crosscheck3.py CERT [trials]"""
import sys, itertools
from fractions import Fraction
import numpy as np
import verify_independent as mv

path = sys.argv[1]
trials = int(sys.argv[2]) if len(sys.argv) > 2 else 6
n, a, t, L, match, fmt, qa, qb = mv.read_cert(path)
assert n == 3
perms, c, idx_id, dcomp, signs, forced = mv.setup_copy(match, a)
QA = np.array(qa, dtype=object).reshape((15,) * n) if fmt == 'full' else mv.expand_orbits(qa, n, perms['tr'])
QB = np.array(qb, dtype=object).reshape((15,) * n) if fmt == 'full' else mv.expand_orbits(qb, n, perms['tr'])
cvec = np.array(c, dtype=object)
Wf = np.multiply.outer(np.multiply.outer(cvec, cvec), cvec).copy()
Wf[(idx_id,) * 3] -= t
import math
lam1 = math.lcm(*[x.denominator for x in Wf.ravel()])
Wint = np.vectorize(lambda x: int(x * lam1), otypes=[object])(Wf)
sL, sR, gA, gB = perms['sL'], perms['sR'], perms['gA'], perms['gB']
QAp = QA + mv.apply_perm(QA, sL) + mv.apply_perm(QA, sR) + mv.apply_perm(mv.apply_perm(QA, sL), sR)
Bint = 4 * L * Wint - lam1 * mv.apply_perm(QAp, gA) - 4 * lam1 * mv.apply_perm(QB, gB)
Zint = Bint + mv.apply_perm(Bint, sL) + mv.apply_perm(Bint, sR) + mv.apply_perm(mv.apply_perm(Bint, sL), sR)
tofl = lambda X, s: np.vectorize(lambda v: float(Fraction(int(v), s)))(X)
QAf, QBf, Zf = tofl(QAp, 4 * L), tofl(QB, L), tofl(Zint, 16 * L * lam1)
D = np.array([mv.diagram_matrix(m) for m in match], dtype=float)


def expect(X, psi):        # psi: per-copy grouped tensor (27, 27, 27)
    U = np.einsum('cmn,ijn->cijm', D, psi)
    V = np.einsum('blj,cijm->bcilm', D, U)
    R = np.einsum('klm,bcilm->bcik', psi.conj(), V)
    M = np.einsum('aki,bcik->abc', D, R)
    return np.sum(X * M)


def group(x, y1, y2):      # x (x) y1 (x) y2 on (A1A2A3)(B1B2B3)(B'1B'2B'3) -> (A1B1B'1)(A2B2B'2)(A3B3B'3)
    T = np.einsum('i,j,k->ijk', x, y1, y2).reshape([3] * 9)
    return T.transpose(0, 3, 6, 1, 4, 7, 2, 5, 8).reshape(27, 27, 27)


Pk = np.zeros((3, 3, 3, 3))   # Pk[k, l, i, j] = <k|P_a(|i><j|)|l>
af = float(a)
for k, l, i, j in itertools.product(range(3), repeat=4):
    Pk[k, l, i, j] = af * (i == j) * (k == l) - ((k == i) * (l == j) + (k == j) * (l == i)) / 4


def map3(rho):             # P_a^{(x)3} applied copy by copy to a 27 x 27 matrix
    T = rho.reshape(3, 3, 3, 3, 3, 3)
    T = np.einsum('klij,ibcjef->kbclef', Pk, T)
    T = np.einsum('klij,aicdjf->akcdlf', Pk, T)
    T = np.einsum('klij,abidej->abkdel', Pk, T)
    return T.reshape(27, 27)


J = np.zeros((729, 729))   # Choi matrix of P_a^{(x)3}, inputs (A1A2A3) then outputs (B1B2B3)
for I in range(27):
    for Jx in range(27):
        E = np.zeros((27, 27)); E[I, Jx] = 1
        J[I * 27:(I + 1) * 27, Jx * 27:(Jx + 1) * 27] = map3(E)

rng = np.random.default_rng(7)
worst = 0.0
cases = []
for _ in range(trials):
    cases.append((rng.normal(size=27) + 1j * rng.normal(size=27), rng.normal(size=27) + 1j * rng.normal(size=27)))
e1 = np.zeros(3); e1[0] = 1
p = np.kron(np.kron(e1, e1), e1).astype(complex)
cases.append((p, p))
v = rng.normal(size=3); v /= np.linalg.norm(v); q = np.kron(np.kron(v, v), v).astype(complex)
cases.append((q, q))
for x, y in cases:
    ny, nx = np.vdot(y, y).real, np.vdot(x, x).real
    val_map = np.vdot(y, map3(np.outer(x.conj(), x.conj().conj())) @ y).real
    val_choi = np.vdot(np.kron(x, y), J @ np.kron(x, y)).real
    lhs = val_map * ny - t.numerator / t.denominator * nx * ny ** 2
    tz = expect(Zf, group(x, y, y)).real
    ta = expect(QAf, group(x.conj(), y, y)).real
    tb = expect(QBf, group(x, y.conj(), y)).real
    rhs = tz + ta + tb
    rel = abs(lhs - rhs) / (nx * ny ** 2)
    worst = max(worst, rel, abs(val_map - val_choi) / (nx * ny))
    print(f'LHS {lhs:+.6e}  RHS {rhs:+.6e}  (Z {tz:+.3e}, QA\' {ta:+.3e}, QB {tb:+.3e})  map vs Choi {val_map - val_choi:+.1e}  rel.diff {rel:.1e}', flush=True)
print(f'largest relative discrepancy {worst:.1e}')
