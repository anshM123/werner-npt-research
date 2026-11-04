"""Fast sector function mu_c = lambda_min of Kt_c / 2^(n-c) built directly from the explicit sector formula
    Kt_c = I + sum_{s vacuum} Abar_s - 1/2 sum_{s excited} |Ut_s><Ut_s| + 1/2 sum_{t excited, s vacuum} Z_s^* Z_t
(single-excitation range, gauge sum A_s = I).  Validated against povm_form.sector_mu."""
import numpy as np, itertools
from povm_form import sqrt_psd

class Sector:
    def __init__(self, n, c):
        self.n, self.c = n, c
        self.strs = [s for s in itertools.product(range(3), repeat=n) if sum(x != 0 for x in s) == c]
        self.pos = {s: i for i, s in enumerate(self.strs)}; self.m = len(self.strs)
        # precompute index lists for the three kinds of terms
        self.vac = [[t for t in range(n) if s[t] == 0] for s in self.strs]
        self.flip = []   # (row string i, col string j, copy t, letter_row, letter_col): same string except the letter at t
        self.hop = []    # (row i (letter z at s), col j (letter z' at t)), s vacuum in col, t vacuum in row
        for j, s in enumerate(self.strs):
            for t in range(n):
                if s[t] == 0:
                    continue
                for zp in (1, 2):
                    sp = list(s); sp[t] = zp; self.flip.append((self.pos[tuple(sp)], j, t, zp, s[t]))
                for u in range(n):
                    if s[u] != 0:
                        continue
                    for z in (1, 2):
                        sp = list(s); sp[t] = 0; sp[u] = z; self.hop.append((self.pos[tuple(sp)], j, u, z, t, s[t]))
    def mu(self, A):
        n, m = self.n, self.m; U = [sqrt_psd(a) for a in A]; K = np.zeros((2 * m, 2 * m), complex)
        for i, vac in enumerate(self.vac):
            B = np.eye(2) + sum((A[t].T for t in vac), np.zeros((2, 2)))      # (I + sum_vac Abar)[r, r'] = delta + sum A[r', r]
            K[i::m, i::m] += B
        for (i, j, t, zrow, zcol) in self.flip:                               # -1/2 U_t[zrow, r] conj(U_t[zcol, r'])
            K[i::m, j::m] += -0.5 * np.outer(U[t][zrow - 1, :], U[t][zcol - 1, :].conj())
        for (i, j, s, z, t, zt) in self.hop:                                  # +1/2 U_s[z, r] conj(U_t[zt, r'])
            K[i::m, j::m] += 0.5 * np.outer(U[s][z - 1, :], U[t][zt - 1, :].conj())
        return np.linalg.eigvalsh((K + K.conj().T) / 2)[0] / 2 ** (n - self.c)
