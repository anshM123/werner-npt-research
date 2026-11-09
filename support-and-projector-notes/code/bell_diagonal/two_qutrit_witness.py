"""Two-qutrit Bell-diagonal witness (Section 7, "Bell-diagonal operators").
K_h = sum_{z in F_3^4} h(z) |b_z><b_z| on C^9 (x) C^9, b_z = vec(D_z)/3.  The weights h in data/bell_k2_witness_h.txt
satisfy all 1440 stabilizer-plane conditions (for every Lagrangian L of F_3^4, every coset L + x and every index-3
subgroup I of L, the three sums of h over the cosets of I in L + x lie in the cone sum >= ||.||_2), yet the rank-two
operator C in data/bell_k2_witness_C_*.txt has <vec C|K_h|vec C> < 0.  This script recomputes both facts from the data.
usage: python two_qutrit_witness.py"""
import itertools, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); data = os.path.join(here, '..', '..', 'data')
h = np.loadtxt(os.path.join(data, 'bell_k2_witness_h.txt'))
C = np.loadtxt(os.path.join(data, 'bell_k2_witness_C_real.txt')) + 1j * np.loadtxt(os.path.join(data, 'bell_k2_witness_C_imag.txt'))
om = np.exp(2j * np.pi / 3)
X1 = np.roll(np.eye(3), 1, axis=0); Z1 = np.diag([1, om, om ** 2])
pts = [np.array(p) for p in itertools.product(range(3), repeat=4)]
key = {tuple(p): i for i, p in enumerate(pts)}
def Dop(p):
    a1, a2, b1, b2 = p
    return np.kron(np.linalg.matrix_power(X1, a1) @ np.linalg.matrix_power(Z1, b1),
                   np.linalg.matrix_power(X1, a2) @ np.linalg.matrix_power(Z1, b2))
Dm = np.array([Dop(p) for p in pts])
def symp(u, v):
    return (u[0] * v[2] + u[1] * v[3] - u[2] * v[0] - u[3] * v[1]) % 3
nz = [p for p in pts if p.any()]
lagr = []
for u in nz:
    for v in nz:
        if symp(u, v): continue
        span = frozenset(tuple((i * u + j * v) % 3) for i in range(3) for j in range(3))
        if len(span) == 9 and span not in lagr: lagr.append(span)
margins = []
for Ls in lagr:
    subs = []
    for w in Ls:
        w = np.array(w)
        if not w.any(): continue
        sp = frozenset(tuple((c * w) % 3) for c in range(3))
        if sp not in subs: subs.append(sp)
    seen = set()
    for x in pts:
        Lx = frozenset(tuple((np.array(z) + x) % 3) for z in Ls)
        if Lx in seen: continue
        seen.add(Lx)
        for I in subs:
            rem = set(Lx); sums = []
            while rem:
                z0 = np.array(next(iter(rem))); g = {tuple((z0 + np.array(i)) % 3) for i in I}
                sums.append(sum(h[key[z]] for z in g)); rem -= g
            s = np.array(sums); margins.append(s.sum() - np.linalg.norm(s))
print(f"Lagrangian planes: {len(lagr)}; stabilizer conditions: {len(margins)}; smallest margin {min(margins):.3e}")
sv = np.linalg.svd(C, compute_uv=False)
v = C.reshape(-1)
K = sum(h[i] * np.outer(Dm[i].reshape(-1), Dm[i].reshape(-1).conj()) / 9 for i in range(81))
print(f"rank of C: {(sv > 1e-9 * sv[0]).sum()} (singular values {np.round(sv[:3] / np.linalg.norm(C), 6)})")
print(f"<vec C|K_h|vec C> / ||C||^2 = {(v.conj() @ K @ v).real / np.vdot(v, v).real:.6e}  (negative: K_h is not 2-block-positive)")
c = np.einsum('zij,ij->z', Dm.conj(), C); p = np.abs(c) ** 2 / (9 * np.vdot(C, C).real)
print(f"Pauli distribution of C: {(p > 1e-3).sum()} points above 1e-3 (not supported on a Lagrangian coset)")
