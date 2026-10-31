"""Section 1: the identity for the aligned class, checked on random members.
   y_k = |k> (x) a + |2> (x) alpha_k,  z_k = |k> (x) b + |2> (x) beta_k  (k = 0, 1),  C = sum_k |y_k><z_k|,  X = sum_k |alpha_k><beta_k|:
     q_n(t, C) = (1-t) q_{n-1}(X) + (2-4t) q_{n-1}(|a><b|) + (1-4t^2) sum_k q_{n-1}(|a><beta_k|)
                 + sum_k <v_k, W v_k>,   v_k = alpha_k (x) b - 2t a (x) beta_k,   W = prod_s (I - t F_s).
   Also checks the crossed-pairing formula and applyB against the defining sum."""
import numpy as np
from werner import d, q, q_pair, q_direct, W_apply

rng = np.random.default_rng(0)


def cv(k):
    return rng.standard_normal(d ** k) + 1j * rng.standard_normal(d ** k)


print("crossed pairing q(|x><y|, |u><v|) = <x (x) v|W|u (x) y>, and applyB against the defining sum (n = 3, t = 0.37):")
x, y, u, v = cv(3), cv(3), cv(3), cv(3)
lhs = q_pair(np.outer(x, y.conj()), np.outer(u, v.conj()), 3, 0.37)
rhs = np.vdot(np.kron(x, v), W_apply(np.kron(u, y), 3, 0.37))
X = rng.standard_normal((27, 27)) + 1j * rng.standard_normal((27, 27))
print(f"  pairing difference {abs(lhs - rhs):.1e};  q(X) via B minus defining sum {abs(q(X, 3, 0.37) - q_direct(X, 3, 0.37)):.1e}")

print("aligned class: largest relative error of the identity, and the smallest q_n/||C||^2 seen (20 random members)")
for n in (2, 3, 4):
    for t in (0.1, 0.3, 0.45, 0.5):
        err, low = 0.0, np.inf
        for trial in range(20):
            a, b = cv(n - 1), cv(n - 1)
            al, be = [cv(n - 1), cv(n - 1)], [cv(n - 1), cv(n - 1)]
            C = np.zeros((d ** n, d ** n), complex)
            for k in range(2):
                ek = np.eye(d)[k]; e2 = np.eye(d)[2]
                C += np.outer(np.kron(ek, a) + np.kron(e2, al[k]), (np.kron(ek, b) + np.kron(e2, be[k])).conj())
            Xc = sum(np.outer(al[k], be[k].conj()) for k in range(2))
            rhs = (1 - t) * q(Xc, n - 1, t) + (2 - 4 * t) * q(np.outer(a, b.conj()), n - 1, t)
            rhs += (1 - 4 * t * t) * sum(q(np.outer(a, be[k].conj()), n - 1, t) for k in range(2))
            for k in range(2):
                vk = np.kron(al[k], b) - 2 * t * np.kron(a, be[k])
                rhs += np.vdot(vk, W_apply(vk, n - 1, t)).real
            nC = np.vdot(C, C).real
            err = max(err, abs(q(C, n, t) - rhs) / nC); low = min(low, q(C, n, t) / nC)
        print(f"  n={n} t={t}: relative error {err:.1e};  min q_n/||C||^2 = {low:+.4f}")
