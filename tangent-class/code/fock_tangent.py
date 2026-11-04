"""Tangent-class operators in Fock variables (F = sum_r |r> (x) D^{1/2} v_r), ordering index = r*3^n + string.
Independent q_direct for checking."""
import numpy as np, itertools
d = 3
def strings(n): return list(itertools.product(range(d), repeat=n))
def q_direct(X, n):
    T = X.reshape([d] * (2 * n)); tot = 0.0
    for J in itertools.product([0, 1], repeat=n):
        Y = T; cur = list(range(n))
        for s in sorted([s for s in range(n) if J[s]], reverse=True):
            i = cur.index(s); m = len(cur); Y = np.trace(Y, axis1=i, axis2=m + i); cur.pop(i)
        tot += (-0.5) ** sum(J) * np.sum(abs(Y) ** 2)
    return tot
def tangent_range(U, a, n):
    st = strings(n); idx = {s: i for i, s in enumerate(st)}; u = np.zeros((2, d ** n), complex)
    for r in range(2):
        u[r, 0] += a[r]
        for s in range(n):
            for z in (1, 2):
                t = [0] * n; t[s] = z; u[r, idx[tuple(t)]] += U[s][z - 1, r]
    return u
class Ops:
    def __init__(self, U, a):
        n = len(U); self.n = n; self.U = [np.asarray(x, complex) for x in U]; self.a = np.asarray(a, complex)
        st = strings(n); self.st = st; idx = {s: i for i, s in enumerate(st)}; D = d ** n; self.D = D
        self.exc = np.array([sum(x != 0 for x in s) for s in st])
        self.Dh = np.array([np.prod([np.sqrt(0.5) if c == 0 else 1.0 for c in s]) for s in st])
        Ab = [(u.conj().T @ u).T for u in self.U]; self.Ab = Ab
        # one-body: sum_u Abar_u (x) (2|0><0|_u + Pi_E_u)
        I2 = np.eye(2)
        self.AbU = np.kron(sum(Ab), np.eye(D))
        Vac = np.zeros((2 * D, 2 * D), complex)
        for i, s in enumerate(st):
            M = sum((Ab[w] for w in range(n) if s[w] == 0), np.zeros((2, 2)))
            for r in range(2):
                for rp in range(2): Vac[r * D + i, rp * D + i] += M[r, rp]
        self.Vac = Vac
        # beta_s : F -> strings (no ref), (beta_s F)(rho) = sum_{z,r} conj(U_s[z,r]) F_r(rho + z@s), rho_s = 0
        B = []
        for s in range(n):
            M = np.zeros((D, 2 * D), complex)
            for j, t in enumerate(st):
                if t[s] == 0: continue
                rho = list(t); rho[s] = 0; i = idx[tuple(rho)]
                for r in range(2): M[i, r * D + j] += np.conj(self.U[s][t[s] - 1, r])
            B.append(M)
        self.B = B; self.T = sum(b.conj().T @ b for b in B); self.al = sum(B)
        self.abra = np.kron(self.a.conj()[None, :], np.eye(D))     # <a| (x) I : F -> strings
        self.Afr = self.al - np.sqrt(2) * self.abra                   # displaced annihilation
        self.Q = self.AbU + self.Vac - self.T + 0.5 * self.Afr.conj().T @ self.Afr
        self.Q = (self.Q + self.Q.conj().T) / 2
    def F_of_v(self, v): return np.concatenate([self.Dh * v[r] for r in range(2)])
    def sector_idx(self, c):
        i = np.where(self.exc == c)[0]; return np.concatenate([i, self.D + i])
def gauge(U, a):
    S = np.outer(np.conj(a), a) + sum(u.conj().T @ u for u in U); w, V = np.linalg.eigh(S); T = V @ np.diag(w ** -0.5) @ V.conj().T
    return [u @ T for u in U], T.T @ a
def rand_tangent(rng, n, scale_a=1.0):
    U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]
    a = (rng.normal(size=2) + 1j * rng.normal(size=2)) * scale_a
    return gauge(U, a)
