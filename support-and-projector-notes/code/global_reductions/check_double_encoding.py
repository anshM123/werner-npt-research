"""Doubly balanced encoding C -> E_r C F_r^*:
q_{k+8r}(E C F^*) = (1/9)[q_k(C) + 2^{-3r} sum_{nu != 0} q_k(C U_nu)], every vector of ran E, ran F locally maximally mixed.
Direct sparse evaluation of q on k + 8r copies (k = 1, 2; r = 1)."""
import itertools
import numpy as np
from qtools import q_dense, q_sparse_np, local_weyl, rand_c, marginals_sparse

LAB = [(a, b) for a in range(3) for b in range(3)]
ELL = [l for l in LAB if l != (0, 0)]


def code(lam, r, shift=0):
    return tuple(((l[0] * lam[0] + l[1] * lam[1] + shift) % 3) for l in ELL for _ in range(r))


def word(i, k):
    return tuple(i // 3 ** (k - 1 - s) % 3 for s in range(k))


def encode_vec(v, k, r, shift):
    out = {}
    for lam in LAB:
        u = local_weyl(k, lam) @ v
        f = code(lam, r, shift)
        for i in range(3 ** k):
            if abs(u[i]) > 0:
                out[word(i, k) + f] = out.get(word(i, k) + f, 0) + u[i] / 3
    return out


rng = np.random.default_rng(11)
ok = True
for k, r in ((1, 1), (2, 1)):
    D = 3 ** k
    x, b, a, y = (rand_c(D, rng) for _ in range(4))
    C = np.outer(x, b.conj()) + np.outer(a, y.conj())
    R, Cc, V = [], [], []
    for lam in LAB:
        for mu in LAB:
            B = local_weyl(k, lam) @ C @ local_weyl(k, mu).conj().T / 9
            fr, fc = code(lam, r, 0), code(mu, r, 1)
            for i in range(D):
                for j in range(D):
                    R.append(word(i, k) + fr); Cc.append(word(j, k) + fc); V.append(B[i, j])
    lhs = q_sparse_np(np.array(R), np.array(Cc), np.array(V))
    tail = sum(q_dense(C @ local_weyl(k, nu), k) for nu in LAB if nu != (0, 0))
    rhs = (q_dense(C, k) + 2.0 ** (-3 * r) * tail) / 9
    dev = abs(lhs - rhs) / max(1, abs(rhs))
    print(f'k={k} r={r} copies={k + 8 * r}: q direct {lhs.real:.12f}  formula {rhs.real:.12f}  rel dev {dev:.1e}')
    ok &= dev < 1e-10
    # marginals of encoded vectors (row code and shifted column code)
    for shift in (0, 1):
        ev = encode_vec(rand_c(D, rng), k, r, shift)
        nrm = sum(abs(t) ** 2 for t in ev.values())
        bad = max(np.abs(m - nrm / 3 * np.eye(3)).max() for m in marginals_sparse(ev))
        ok &= bad < 1e-10
        print(f'   shift {shift}: max marginal deviation from I/3: {bad:.1e}')
# the 81 x 81 flag compression is diagonal with weights 1 (lambda = mu) and 2^{-3r}, r = 1, 2, 3
for r in (1, 2, 3):
    rows = [code(l, r, 0) for l in LAB]
    cols = [code(m, r, 1) for m in LAB]
    states = [(rows[i], cols[j], i == j) for i in range(9) for j in range(9)]
    bad = 0
    for (ra, ca, same) in states:
        for (rb, cb, _) in states:
            val = 1.0
            for s in range(8 * r):
                val *= (ra[s] == rb[s] and ca[s] == cb[s]) - 0.5 * (ra[s] == ca[s] and rb[s] == cb[s])
                if val == 0:
                    break
            exp = (1.0 if same else 2.0 ** (-3 * r)) if (ra == rb and ca == cb) else 0.0
            bad += abs(val - exp) > 0
    print(f'flag compression r={r}: entries off the stated pattern: {bad}')
    ok &= bad == 0
print('PASS' if ok else 'FAIL')
