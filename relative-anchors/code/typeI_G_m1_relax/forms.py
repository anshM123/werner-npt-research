"""Explicit formulas (generic number type via numpy float here) for the type I block of (R_inf)_{k-1}, m = 1, gauge
x_t = 0, after the exact elimination of q:
  Phi_q(X, rho, zeta) = F(X) + <Hfull(rho, zeta), X> + Acal(rho, zeta)
with F the (C_u) form of Theorem C (q eliminated), Acal = 2 rho^2 + |zeta|^2 - m^2/2,
m = kappa rho + w_t (w.zeta)/(sig + S), kappa = (2 sig + S)/(sig + S),
Hfull[r,h] = (rho - m) d_rh + m w_r w_h/(2K) - c_r W_tr w_h + w_t c_r W_hr - (Wy)_t w_r w_h/2 + w_t w_r (Wy)_h/2,
c_r = Y0/2 + w_r zeta_r, Y0 = w_t rho + w.zeta, W = E_0^{-1} on U+t, y = rho e_t + zeta."""
import numpy as np


class Forms:
    def __init__(self, st, sU):
        self.st = float(st); self.sU = np.array(sU, float); u = self.u = len(sU)
        self.w = np.sqrt(self.sU); self.wt = np.sqrt(self.st)
        S = self.S = self.sU.sum(); sig = self.st
        self.K = sig + S / 2
        self.kappa = (2 * sig + S) / (sig + S)
        # E (pattern T, on U) and E0 (pattern empty, on U+t; index 0 = t)
        self.d = self.sU + S / 2 + sig
        self.E = np.diag(self.d) + np.outer(self.w, self.w) / 2
        self.Ei = np.linalg.inv(self.E)
        wtil = np.concatenate([[self.wt], self.w]); self.wtil = wtil
        d0 = np.concatenate([[sig], self.sU]) + (sig + S) / 2; self.d0 = d0
        self.E0 = np.diag(d0) + np.outer(wtil, wtil) / 2
        self.W = np.linalg.inv(self.E0)

    def Hfull(self, rho, zeta):
        w, wt, W = self.w, self.wt, self.W
        Y0 = wt * rho + w @ zeta
        m = self.kappa * rho + wt * (w @ zeta) / (self.st + self.S)
        c = Y0 / 2 + w * zeta
        y = np.concatenate([[rho], zeta])
        Wy = W @ y
        Wtr = W[0, 1:]; WUU = W[1:, 1:]
        H = (rho - m) * np.eye(self.u) + m * np.outer(w, w) / (2 * self.K)
        H += -np.outer(c * Wtr, w) + wt * (c[:, None] * WUU.T)
        H += -Wy[0] * np.outer(w, w) / 2 + wt * np.outer(w, Wy[1:]) / 2
        return H

    def Acal(self, rho, zeta):
        m = self.kappa * rho + self.wt * (self.w @ zeta) / (self.st + self.S)
        return 2 * rho ** 2 + zeta @ zeta - m ** 2 / 2

    def Gcal(self, a, c):
        """Theorem C's representative: F(X'+Omega) - F(X') - F(Omega) = sig <Gcal, X'> for hard X'."""
        Ei, w = self.Ei, self.w
        p = w @ c
        return 2 * Ei @ np.diag(a * w) + 2 * np.diag(c * w) @ Ei + p * (Ei - np.eye(self.u) / self.K)

    def f(self, a, c):
        Ei, w, sig = self.Ei, self.w, self.st
        E0 = self.E - sig * np.eye(self.u)
        z = Ei @ w
        p = w @ c
        return a @ E0 @ Ei @ a + c @ E0 @ Ei @ c - sig * p ** 2 / (2 * self.K ** 2) + p * (a @ z) + 2 * np.sum(a * c * w * z)
