"""Section 9: two copies, positive semidefinite D of any rank and any local dimensions:
4 q_2(D) - (2 Tr D^2 - (Tr D)^2) = 8 Tr[(D (x) D)(P_A (x) P_A)] >= 0.  Check the identity and the sign on random D."""
import numpy as np
rng = np.random.default_rng(5)
def q2(C, d1, d2):
    T = C.reshape(d1, d2, d1, d2)
    t1 = np.einsum('ajak->jk', T); t2 = np.einsum('iaja->ij', T); tr = np.trace(C)
    return np.linalg.norm(C)**2 - 0.5*np.linalg.norm(t1)**2 - 0.5*np.linalg.norm(t2)**2 + 0.25*abs(tr)**2
worst = 1e9; maxdev = 0
for trial in range(400):
    d1, d2 = rng.integers(2, 5, size=2); D = d1*d2; r = rng.integers(1, D+1)
    G = rng.standard_normal((D, r)) + 1j*rng.standard_normal((D, r)); Dm = G @ G.conj().T; Dm /= np.trace(Dm).real
    lhs = 4*q2(Dm, d1, d2) - (2*np.linalg.norm(Dm)**2 - np.trace(Dm).real**2)
    # P_A (x) P_A on (site1 copy1, site1 copy2) (x) (site2 copy1, site2 copy2): build via swaps
    # Tr[(D x D) P_A1 P_A2], P_A = (I - F)/2 per site: expand
    t_I = np.trace(Dm).real**2
    t_F1 = np.trace(np.einsum('ajbj->ab', Dm.reshape(d1, d2, d1, d2)) @ np.einsum('ajbj->ab', Dm.reshape(d1, d2, d1, d2))).real
    t_F2 = np.trace(np.einsum('iaib->ab', Dm.reshape(d1, d2, d1, d2)) @ np.einsum('iaib->ab', Dm.reshape(d1, d2, d1, d2))).real
    t_F12 = np.trace(Dm @ Dm).real
    w12 = (t_I - t_F1 - t_F2 + t_F12) / 4
    maxdev = max(maxdev, abs(lhs - 8*w12)); worst = min(worst, lhs)
print(f"max |identity residual| = {maxdev:.2e};  smallest 4q_2 - (2TrD^2 - (TrD)^2) over 400 random PSD D (ranks 1..D, dims 2..4) = {worst:.3e}")
