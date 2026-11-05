"""The (C_u) form.  Principal block of the type I block (units of 2Y) on the states P_{r,h} (r, h in U) and Q, with q
eliminated.  Parameters: the U-weights w (length u) and sigma = s_T > 0.  Generic in the number type.

Notation: S = sum s_h, d_h = s_h + S/2 + sigma, D = diag(d), E = D + w w^T/2, E_0 = E - sigma I, K = sigma + S/2,
eta = w/d, g = sum s_h/d_h, rho_h = s_h/d_h, z = E^{-1} w = eta/(1 + g/2).

Closed form (derived in NOTES.md, Section 2):
  F(X) = 1/2 ||X||^2 + sum_h X_hh^2 + 1/2 (tr X)^2 - 1/2 (tr X - w^T X w/(2K))^2 - B(X),
  B(X) = sum_{r,r'} < T_{r'} psi_r , E^{-1} T_r psi_{r'} >,  T_{r'} psi = 1/2 (w.psi) e_{r'} + w_{r'} (psi/2 + psi_{r'} e_{r'}).
"""
from fractions import Fraction as Fr


def inv_mat(M, F):
    n = len(M)
    A = [list(M[i]) + [F(1) if i == j else F(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        r = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[r] = A[r], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


class CU:
    def __init__(self, w, sigma, F=Fr):
        self.F = F
        self.w = [F(x) for x in w]
        self.u = len(w)
        self.s = [x * x for x in self.w]
        self.sig = F(sigma)
        u = self.u
        self.S = sum(self.s, F(0))
        self.K = self.sig + self.S / 2
        self.d = [self.s[h] + self.S / 2 + self.sig for h in range(u)]
        self.E = [[(self.d[i] if i == j else F(0)) + self.w[i] * self.w[j] / 2 for j in range(u)] for i in range(u)]
        self.Ei = inv_mat(self.E, F)

    def T1(self, rp, psi):
        """T_{(r',1)} psi on pattern T (free set U)."""
        F = self.F
        wpsi = sum((a * b for a, b in zip(self.w, psi)), F(0))
        out = [self.w[rp] * x / 2 for x in psi]
        out[rp] += wpsi / 2 + self.w[rp] * psi[rp]
        return out

    def form(self, X):
        """F(X) (q eliminated), X a u x u list of lists."""
        F = self.F
        u = self.u
        w = self.w
        A = sum((X[r][h] ** 2 for r in range(u) for h in range(u)), F(0)) / 2
        A += sum((X[h][h] ** 2 for h in range(u)), F(0))
        tr = sum((X[h][h] for h in range(u)), F(0))
        A += tr * tr / 2
        beta = sum((w[r] * X[r][h] * w[h] for r in range(u) for h in range(u)), F(0))
        L = tr - beta / (2 * self.K)
        A -= L * L / 2
        Bv = F(0)
        V = [[self.T1(rp, X[r]) for rp in range(u)] for r in range(u)]   # V[r][rp] = T_{rp} psi_r
        for r in range(u):
            for rp in range(u):
                a = V[r][rp]
                b = V[rp][r]
                Bv += sum((a[i] * self.Ei[i][j] * b[j] for i in range(u) for j in range(u)), F(0))
        return A - Bv

    def full_matrix(self):
        """The (u^2 + 1) x (u^2 + 1) matrix of the principal block (before eliminating q), from the closed form
        with q: Phi = 1/2||X||^2 + sum X_hh^2 + 1/2 tr^2 + 2 q (-tr X + w^T X w/(2K)) + 2 q^2 - B(X)."""
        F = self.F
        u = self.u
        n = u * u + 1
        M = [[F(0)] * n for _ in range(n)]

        def quad(X):
            Fx = self.form(X)
            return Fx

        # polarisation of the q-eliminated form is not what we want here; build Phi directly
        def phi(X, q):
            tr = sum((X[h][h] for h in range(u)), F(0))
            beta = sum((self.w[r] * X[r][h] * self.w[h] for r in range(u) for h in range(u)), F(0))
            L = -tr + beta / (2 * self.K)
            val = self.form(X) + (tr - beta / (2 * self.K)) ** 2 / 2   # = A-part - B
            return val + 2 * q * L + 2 * q * q

        basis = []
        for r in range(u):
            for h in range(u):
                X = [[F(0)] * u for _ in range(u)]
                X[r][h] = F(1)
                basis.append((X, F(0)))
        basis.append(([[F(0)] * u for _ in range(u)], F(1)))

        def add(x, y):
            return ([[a + b for a, b in zip(ra, rb)] for ra, rb in zip(x[0], y[0])], x[1] + y[1])
        diag = [phi(*b) for b in basis]
        for i in range(n):
            M[i][i] = diag[i]
            for j in range(i + 1, n):
                v = phi(*add(basis[i], basis[j]))
                M[i][j] = M[j][i] = (v - diag[i] - diag[j]) / 2
        return M
