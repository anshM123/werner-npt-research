"""27-branch balanced projector: Q = H Pi H^*, N = k + 8r + 2 (r even).
q_N(Q) = (2/27) q_k(Pi) + ((2 eps + eta)/27) S(Pi) + (zeta/81) q_k(T(Pi)).
Exact check of the 729 x 729 flag compression (r = 2), contraction with core data (k = 1, 2), marginals."""
import numpy as np
from fractions import Fraction as Fr
from qtools import q_dense, local_weyl, rand_c, marginals_sparse, OMEGA

LAB = [(a, b) for a in range(3) for b in range(3)]
ELL = [l for l in LAB if l != (0, 0)]
r = 2
BR = [(lam, s) for lam in LAB for s in range(3)]          # 27 branch words


def code(lam, s):
    return tuple(((l[0] * lam[0] + l[1] * lam[1] + s) % 3) for l in ELL for _ in range(r)) + (s, s)


CODES = [code(lam, s) for (lam, s) in BR]
L = len(CODES[0])
states = [(i, j) for i in range(27) for j in range(27)]
eps, eta, zeta = Fr(1, 2 ** (3 * r)), Fr(1, 2 ** (2 * r + 2)), Fr(1, 2 ** (8 * r + 2))

# exact flag compression: entry = prod_p (d(row,row')d(col,col') - 1/2 d(row,col) d(row',col'))
A = np.array([CODES[i] for (i, j) in states])
B = np.array([CODES[j] for (i, j) in states])
t2 = (A == B)                                              # row == col, per state and position
F = {}
bad = 0
for u in range(729):
    t1 = (A[u] == A) & (B[u] == B)                         # (729, L)
    tt = t2[u][None, :] & t2                               # row=col in both
    zero = ~(t1 | tt)
    nz = ~zero.any(axis=1)
    for v in np.nonzero(nz)[0]:
        neg = int(np.sum(~t1[v] & tt[v]))                  # factor -1/2
        half = int(np.sum(t1[v] & tt[v]))                  # factor 1/2
        F[(u, v)] = Fr((-1) ** neg, 2 ** (neg + half))
# expected pattern
for (u, v), val in F.items():
    (i, j), (i2, j2) = states[u], states[v]
    (lam, s), (mu, t) = BR[i], BR[j]
    if i == j:
        exp = zeta if i2 == j2 else None
    elif u == v:
        exp = Fr(1) if (s != t and lam == mu) else (eps if s != t else eta)
    else:
        exp = None
    bad += (exp is None) or (val != exp)
nd = sum(1 for (u, v) in F if states[u][0] == states[u][1] and states[v][0] == states[v][1])
print('nonzero flag entries', len(F), ' off-pattern', bad, ' diagonal-pair block entries', nd, '(expect 729)')
ok = bad == 0 and nd == 729

rng = np.random.default_rng(5)
for k in (1, 2):
    D = 3 ** k
    Qm, _ = np.linalg.qr(rand_c((D, 2), rng))
    Pi = Qm @ Qm.conj().T
    U = {lam: local_weyl(k, lam) for lam in LAB}
    blocks = {(i, j): U[BR[i][0]] @ Pi @ U[BR[j][0]].conj().T for i in range(27) for j in range(27)}

    def ip(X, Y):   # <vec X, M_k vec Y> by polarisation of q_dense (real symmetric M)
        return (q_dense(X + Y, k) - q_dense(X - Y, k) - 1j * q_dense(X + 1j * Y, k) + 1j * q_dense(X - 1j * Y, k)) / 4

    direct = 0
    for (u, v), val in F.items():
        direct += float(val) * ip(blocks[states[u]], blocks[states[v]])
    direct /= 27 ** 2
    S = sum(q_dense(Pi @ U[nu], k) for nu in LAB if nu != (0, 0))
    T = sum(U[lam] @ Pi @ U[lam].conj().T for lam in LAB)
    formula = (2 / 27) * q_dense(Pi, k) + float(2 * eps + eta) / 27 * S + float(zeta) / 81 * q_dense(T, k)
    dev = abs(direct - formula)
    print(f'k={k}: contraction {direct.real:.15f} formula {formula.real:.15f} dev {dev:.1e}')
    ok &= dev < 1e-12
    # approximation bound |q_N(Q) - (2/27) q_k(Pi)| <= 2^{-(2r+1)}
    ok &= abs(formula - 2 / 27 * q_dense(Pi, k)) <= 2.0 ** (-(2 * r + 1))
    # encoded vectors: isometry and marginals
    vecs = []
    for c in range(2):
        v = Qm[:, c]
        out = {}
        for i, (lam, s) in enumerate(BR):
            w = U[lam] @ v
            for x in range(D):
                key = tuple(x // 3 ** (k - 1 - q) % 3 for q in range(k)) + CODES[i]
                out[key] = out.get(key, 0) + w[x] / (3 * np.sqrt(3))
        vecs.append(out)
    keys = set(vecs[0]) | set(vecs[1])
    G = np.array([[sum(np.conj(a.get(z, 0)) * b.get(z, 0) for z in keys) for b in vecs] for a in vecs])
    ok &= np.allclose(G, np.eye(2))
    sup = {z: (vecs[0].get(z, 0) + OMEGA * vecs[1].get(z, 0)) / np.sqrt(2) for z in keys}
    badm = max(np.abs(m - np.eye(3) / 3).max() for m in marginals_sparse(sup))
    print(f'   isometry ok: {np.allclose(G, np.eye(2))}; superposition marginal deviation from I/3: {badm:.1e}')
    ok &= badm < 1e-12
print('PASS' if ok else 'FAIL')
