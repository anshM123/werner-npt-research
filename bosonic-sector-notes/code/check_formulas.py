"""Section 1 checks.
   (a) the annihilation formula for q_n against the partial-trace definition (n = 2, 3, 4);
   (b) the SU(3) spectrum: the form has eigenvalue gamma_{n,k} with multiplicity (k+1)^3 (n = 1..6);
   (c) the closed forms gamma_{n,n-1} = -1/(2n), gamma_{n,n-2} = 1/(2(n-1)), gamma_{n,0} = 2^{-n-2}(1 + (-1)^n (2n+3))."""
import numpy as np
from fractions import Fraction
from bosonic import basis, q_and_grad, gamma, sym_embed, q_direct

rng = np.random.default_rng(1)
print("(a) annihilation formula vs partial traces, random rank-2 C, t = 1/2 and t = 0.37")
for n in (2, 3, 4):
    E = sym_embed(n); D = E.shape[1]
    C = (rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2))) @ \
        (rng.standard_normal((2, D)) + 1j * rng.standard_normal((2, D)))
    for t in (0.5, 0.37):
        a, b = q_and_grad(C, n, t)[0], q_direct(E @ C @ E.T, n, t)
        print(f"  n={n} t={t}: formula {a:.10f}  direct {b:.10f}  diff {abs(a - b):.1e}")

print("(b) spectrum of the q-form on End(Sym^n), t = 1/2")
for n in range(1, 7):
    D = len(basis(n)); Q = np.zeros((D * D, D * D))
    for e in range(D * D):
        X = np.zeros(D * D); X[e] = 1
        Q[:, e] = q_and_grad(X.reshape(D, D), n)[1].ravel()
    ev = np.linalg.eigvalsh((Q + Q.T) / 2)
    predicted = sorted(float(gamma(n, k)) for k in range(n + 1) for _ in range((k + 1) ** 3))
    print(f"  n={n}: max |eig - predicted| = {np.abs(np.sort(ev) - np.array(predicted)).max():.1e};  "
          f"gamma_(n,k), k = 0..n: {[str(gamma(n, k)) for k in range(n + 1)]}")

print("(c) closed forms")
ok = True
for n in range(2, 40):
    ok &= gamma(n, n - 1) == Fraction(-1, 2 * n)
    ok &= gamma(n, n - 2) == Fraction(1, 2 * (n - 1))
    ok &= gamma(n, 0) == Fraction(1 + (-1) ** n * (2 * n + 3), 2 ** (n + 2))
    ok &= all((gamma(n, k) > 0) == ((n - k) % 2 == 0) for k in range(n + 1))
print("  closed forms and sign pattern (-1)^(n-k) hold for n = 2..39:", ok)
