"""Float sanity check of the whole argument directly from u1form (no claim rests on it).
For random and extreme weight vectors (m up to 64): Ghat = diag(1, c) G diag(1, c) blockwise, c^2 = 1 + sigma/S_T;
reports  max ||Ghat_tt'||_F S_T / sqrt(s_t s_t')  (proved <= 2/5),
         min lambda_min(Ghat_tt) / (1 - s_t/S_T)  (proved > 2/5),
         max_t sum_{t'} ||Ghat_tt'||_2 sqrt(s_t'/s_t) / lambda_min(Ghat_tt)  (proved < 1),  and min eig of G."""
import numpy as np
from u1form import phi_matrix

rng = np.random.default_rng(11)
worstA = 0.0; worstB = 9.0; worstR = 0.0; mineig = 9.0
for trial in range(400):
    m = int(rng.choice([1, 2, 3, 5, 8, 16, 32, 64]))
    kind = trial % 5
    if kind == 0: w = rng.random(m + 1) + 0.05
    elif kind == 1: w = np.exp(rng.uniform(-6, 2, m + 1))
    elif kind == 2: w = np.ones(m + 1); w[1] = np.sqrt(max(m - 1, 1.0)); w[0] = np.exp(rng.uniform(-4, 4))
    elif kind == 3: w = np.exp(rng.uniform(-4, 0, m + 1)); w[1] = 1.0; w[0] = np.exp(rng.uniform(-5, 5))
    else: w = np.exp(rng.uniform(-1, 1, m + 1)); w[0] = np.exp(rng.uniform(-6, 6))
    P = np.array(phi_matrix(list(w), float)); P = (P + P.T) / 2
    iq = 2 * m
    G = P[:iq, :iq] - np.outer(P[:iq, iq], P[iq, :iq]) / P[iq, iq]
    s = w[1:] ** 2; ST = s.sum(); sig = w[0] ** 2; c = np.sqrt(1 + sig / ST)
    Dv = np.ones(iq); Dv[m:] = c
    Gh = Dv[:, None] * G * Dv[None, :]
    blk = lambda t, u: Gh[np.ix_([t, m + t], [u, m + u])]
    mineig = min(mineig, np.linalg.eigvalsh(G)[0])
    for t in range(m):
        lm = np.linalg.eigvalsh(blk(t, t))[0]
        worstB = min(worstB, lm / (1 - s[t] / ST) if m > 1 else 9.0)
        off = 0.0
        for u in range(m):
            if u == t: continue
            worstA = max(worstA, np.linalg.norm(blk(t, u)) * ST / np.sqrt(s[t] * s[u]))
            off += np.linalg.norm(blk(t, u), 2) * np.sqrt(s[u] / s[t])
        worstR = max(worstR, off / lm)
print(f"max ||Ghat_tt'||_F S_T/sqrt(s_t s_t') = {worstA:.4f} (proved <= 0.4)")
print(f"min lambda_min(Ghat_tt)/(1 - s_t/S_T)  = {worstB:.4f} (proved > 0.4)")
print(f"max weighted dominance ratio           = {worstR:.4f} (proved < 1)")
print(f"min eigenvalue of G                    = {mineig:.3e} (proved > 0)")
print("PASS" if worstA <= 0.4 and worstB > 0.4 and worstR < 1 and mineig > 0 else "FAIL")
