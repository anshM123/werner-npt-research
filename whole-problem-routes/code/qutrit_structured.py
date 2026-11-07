"""Check for Observation 7.4 (whole-problem-routes/README.md): structured qutrit inputs with full one-copy support (qutrit GHZ, eps, AME(4,3), Fourier-basis GHZ real
and imaginary parts, and products of these with the plane witness) for P_a^{(x)n}: least output eigenvalue at
a = a_n + 0.001 and a = 3/4, then alternating descent from each at a = a_n + 0.001 (n <= 5).
usage: python qutrit_structured.py   (log: logs/qutrit_structured.log)"""
import sys, itertools
import numpy as np
def phi_for(a, n):
    d = 3; e = np.eye(d); M = np.zeros((9, 9))
    for i in range(3):
        for j in range(3):
            E = np.outer(e[i], e[j]); M[:, 3 * i + j] = (a * np.trace(E) * np.eye(3) - (E + E.T) / 4).reshape(9)
    M4 = M.reshape(3, 3, 3, 3); D = 3 ** n
    def Phi(X):
        T = X.reshape([3] * (2 * n))
        for c in range(n):
            T = np.moveaxis(np.tensordot(M4, T, axes=([2, 3], [c, n + c])), [0, 1], [c, n + c])
        return T.reshape(D, D)
    return Phi
def ket(digits):
    v = np.zeros(3 ** len(digits)); v[int(''.join(map(str, digits)), 3)] = 1; return v
def ghz3(n): return sum(ket([k] * n) for k in range(3)) / np.sqrt(3)
def ghzy(n):
    p = np.array([1, 1j, 0]) / np.sqrt(2); g = np.ones(1, dtype=complex)
    for _ in range(n): g = np.kron(g, p)
    v = np.real(g + np.conj(g)); return v / np.linalg.norm(v)
def fourier_ghz(n):
    w = np.exp(2j * np.pi / 3); F = [np.array([1, w ** k, w ** (2 * k)]) / np.sqrt(3) for k in range(3)]
    g = 0
    for k in range(3):
        t = np.ones(1, dtype=complex)
        for _ in range(n): t = np.kron(t, F[k])
        g = g + t
    return g / np.linalg.norm(g)
def eps3():
    v = np.zeros(27)
    for p in itertools.permutations(range(3)): v[9 * p[0] + 3 * p[1] + p[2]] = np.linalg.det(np.eye(3)[list(p)])
    return v / np.linalg.norm(v)
def ame4():
    v = np.zeros(81)
    for i in range(3):
        for j in range(3): v[int(''.join(map(str, [i, j, (i + j) % 3, (i + 2 * j) % 3])), 3)] = 1
    return v / 3
def outputs(x, a, n):
    Phi = phi_for(a, n)
    X = np.real(np.outer(x, np.conj(x)))     # the map sees Re(xx*) only
    return np.linalg.eigvalsh(Phi(X))[0]
def descend(x, a, n, iters=500):
    Phi = phi_for(a, n); x = np.real(x) if np.isrealobj(x) else np.real(x)
    x = x / np.linalg.norm(x); val = np.inf
    for it in range(iters):
        w, V = np.linalg.eigh(Phi(np.outer(x, x))); y = V[:, 0]
        w2, V2 = np.linalg.eigh(Phi(np.outer(y, y))); nv = w2[0]; x = V2[:, 0]
        if abs(nv - val) < 1e-15: val = nv; break
        val = nv
    return val
cands = {3: {'ghz3': ghz3(3), 'eps': eps3(), 'fourierGHZ': fourier_ghz(3), 'ReFourierGHZ': np.real(fourier_ghz(3)),
             'ImFourierGHZ': np.imag(fourier_ghz(3))},
         4: {'ghz3': ghz3(4), 'ame4': ame4(), 'fourierGHZ': fourier_ghz(4), 'ReFourierGHZ': np.real(fourier_ghz(4)),
             'ghzy2 x ghz3_2': np.kron(ghzy(2), ghz3(2))},
         5: {'ghz3': ghz3(5), 'eps x ghzy2': np.kron(eps3(), ghzy(2)), 'ReFourierGHZ': np.real(fourier_ghz(5)),
             'ame4 x e0': np.kron(ame4(), np.eye(3)[0])},
         6: {'ghz3': ghz3(6), 'eps x eps': np.kron(eps3(), eps3()), 'ame4 x ghzy2': np.kron(ame4(), ghzy(2)),
             'ReFourierGHZ': np.real(fourier_ghz(6))}}
for n, cs in cands.items():
    an = (1 + 2 ** (1 - 1 / n)) / 4
    for name, x in cs.items():
        x = x / np.linalg.norm(x)
        l1 = outputs(x, an + 0.001, n); l2 = outputs(x, 0.75, n)
        dv = descend(x, an + 0.001, n) if n <= 5 else float('nan')
        print(f'n = {n}, {name:16s}: lambda_min at a_n + 0.001: {l1:+.4e}; at 3/4: {l2:+.4e}; after descent at a_n + 0.001: {dv:+.4e}  (m^n there {(an + 0.001 - 0.5) ** n:.3e})', flush=True)
