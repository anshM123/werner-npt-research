"""Exploration (floats): the chain C1.
  F(X) = F(X') + sigma <G, X'> + sigma f(a, c),  F(X') >= (5/12) sum_{r!=h} X'_rh^2 + sum_h Lam_h X'_hh^2, Lam_h = 3/2 - rho_h,
  so F(X) >= sigma [ f - (sigma/4) ( (12/5) sum_{r!=h} G_rh^2 + sum_h G_hh^2 / Lam_h ) ] =: sigma Phi1(a, c),
with G = 2 E^{-1} diag(a w) + 2 diag(c w) E^{-1} + (w.c)(E^{-1} - I/K) (unprojected representative).
Report min generalized eigenvalue of Phi1 w.r.t. f on {a perp w}; also variants."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def null_basis(A, tol=1e-12):
    u_, s_, vt = np.linalg.svd(A)
    r = (s_ > tol * max(1, s_.max())).sum()
    return vt[r:].T


def mats(w, sig):
    w = np.asarray(w, float); u = len(w); s = w ** 2; S = s.sum(); K = sig + S / 2
    d = s + S / 2 + sig
    E = np.diag(d) + np.outer(w, w) / 2
    Ei = np.linalg.inv(E)
    E0 = E - sig * np.eye(u)
    rho = s / d; g = rho.sum()
    z = Ei @ w
    Bw = null_basis(w[None, :])
    # coordinates y = (alpha (u-1), c (u)), a = Bw alpha
    n = 2 * u - 1
    def ac(y):
        return Bw @ y[:u - 1], y[u - 1:]
    def f(y):
        a, c = ac(y)
        cw = c @ w
        return (a @ E0 @ Ei @ a + c @ E0 @ Ei @ c - sig * cw ** 2 / (2 * K * K) + cw * (a @ z) + 2 * np.sum(a * c * w * z))
    def G(y):
        a, c = ac(y)
        return 2 * Ei @ np.diag(a * w) + 2 * np.diag(c * w) @ Ei + (c @ w) * (Ei - np.eye(u) / K)
    Lam = 1.5 - rho
    def coup(y, mode):
        Gm = G(y)
        if mode == 'proj':
            P = np.eye(u) - np.outer(w, w) / S
            Gm = P @ Gm @ P
        off = (Gm ** 2).sum() - (np.diag(Gm) ** 2).sum()
        return (sig / 4) * ((12 / 5) * off + (np.diag(Gm) ** 2 / Lam).sum())
    def qmat(fun):
        I = np.eye(n)
        M = np.zeros((n, n))
        dg = [fun(I[i]) for i in range(n)]
        for i in range(n):
            M[i, i] = dg[i]
            for j in range(i + 1, n):
                M[i, j] = M[j, i] = (fun(I[i] + I[j]) - dg[i] - dg[j]) / 2
        return M
    Fm = qmat(f)
    C1 = qmat(lambda y: coup(y, 'raw'))
    C2 = qmat(lambda y: coup(y, 'proj'))
    L = np.linalg.cholesky(Fm); Li = np.linalg.inv(L)
    return np.linalg.eigvalsh(Li @ (Fm - C1) @ Li.T)[0], np.linalg.eigvalsh(Li @ (Fm - C2) @ Li.T)[0]


if __name__ == "__main__":
  rng = np.random.default_rng(4)
  for u in (2, 3, 4, 5, 6, 8):
      worst = [np.inf, np.inf]
      arg = [None, None]
      for trial in range(400):
          spread = rng.choice([0.3, 1, 2, 4])
          w = np.exp(rng.uniform(-spread, spread, u))
          if trial % 4 == 1:
              w[0] *= np.exp(rng.uniform(2, 5))
          sig = np.exp(rng.uniform(-7, 7)) * (w ** 2).sum()
          r = mats(w, sig)
          for i in range(2):
              if r[i] < worst[i]:
                  worst[i] = r[i]; arg[i] = (np.round(w ** 2 / (w ** 2).sum(), 3), sig / (w ** 2).sum())
      print(u, 'C1 raw: %.4f   C1 projected: %.4f' % tuple(worst), 'argmin raw', arg[0], flush=True)
