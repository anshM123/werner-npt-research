"""Checks for the local-overlap budget.
(1) Symbolic: beta(alpha, tau) = 1 - v/u for the compressed family X(x) = P(x)Qbar - x|omega><omega| with base x = -1
    (compressed P^s), core x = alpha_0, Gurvits-Barnum separable point x_s = 1/tau, equals
    (1 + alpha_0)(alpha tau - 1)/((alpha_0 tau - 1)(1 + alpha)), and at tau = d equals the Werner budget.
(2) Numerical: for random rank-p, rank-q projections P, Q in C^d with tau = Tr(P Q) = ||(P (x) Qbar) Omega||^2, the compressed operator
    1 - x|omega><omega| on supp is separable-certified (PSD with ||X - 1||_2 <= 1) for |x| tau <= 1.
"""
import numpy as np
import sympy as sp

x, x0, tau, d = sp.symbols('x x0 tau d', positive=True)
# X(x) = lam * base + mu * core in the |omega><omega| coordinate (P (x) Qbar coefficient normalised to 1)
lam = (x0 - x) / (x0 + 1); mu = (1 + x) / (1 + x0)
assert sp.simplify(lam + mu - 1) == 0 and sp.simplify(lam * 1 + mu * (-x0) + x) == 0
u = mu / lam
v = u.subs(x, 1 / tau)
beta = sp.simplify(1 - v / u)
target = (1 + x0) * (x * tau - 1) / ((x0 * tau - 1) * (1 + x))
print('beta - target simplifies to', sp.simplify(beta - target))
kap = lambda a: (1 - a) / (1 + a)
werner = (((d - 1) / (d + 1)) - kap(x)) / (((d - 1) / (d + 1)) - kap(x0))
print('beta(tau=d) - Werner budget simplifies to', sp.simplify(target.subs(tau, d) - werner))
# examples
f = sp.lambdify((x, x0, tau), target)
print('d=3 endpoint core, alpha=5/11, tau=3:', f(5 / 11, 0.5, 3), '(Werner 3/4)')
print('alpha=0.4, tau=2.4 (free since alpha tau <= 1):', max(0.0, f(0.4, 0.5, 2.4)))
print('alpha=0.4, tau=2.8:', f(0.4, 0.5, 2.8), ' tau=3:', f(0.4, 0.5, 3))

rng = np.random.default_rng(3)
worst = 1e9
for trial in range(300):
    D = rng.integers(3, 6)
    p, q = rng.integers(1, D + 1), rng.integers(1, D + 1)
    def proj(r):
        Z = rng.normal(size=(D, r)) + 1j * rng.normal(size=(D, r))
        Qm, _ = np.linalg.qr(Z)
        return Qm
    VP, VQ = proj(p), proj(q)
    om = np.zeros(D * D, dtype=complex)
    for k in range(D):
        e = np.zeros(D); e[k] = 1
        om += np.kron(e, e)
    # omega = (P (x) Qbar) Omega restricted to supp: coordinates in basis VP (x) conj(VQ)
    W = np.kron(VP, VQ.conj())
    w = W.conj().T @ om
    t = np.vdot(w, w).real
    PQ = (VP @ VP.conj().T) @ (VQ @ VQ.conj().T)        # ||(P (x) Qbar) Omega||^2 = Tr(P Qbar^T) = Tr(P Q)
    for xx in np.linspace(-1 / t, 1 / t, 7):        # GB ball; the base x = -1 is separable as a compression of P^s
        X = np.eye(p * q) - xx * np.outer(w, w.conj())
        worst = min(worst, np.linalg.eigvalsh(X)[0])
        assert abs(xx) * t <= 1 + 1e-12
    assert abs(t - np.trace(PQ).real) < 1e-9, (t, np.trace(PQ))
print('random compressions: ||omega||^2 = Tr(P Q) confirmed; min eigenvalue on |x| tau <= 1:', worst)
