"""Section 3: the stored rank-two operator C on (C^3)^{(x)3} (data/sci_counterexample_n3.npy) against the contraction
   inequality R_3(C) >= (1/2) min_{s, w, w'} R_2(<w|_s C|w'>_s) at t = 1/2, with R_k = q_k/||.||^2.
   The inner minimum is computed by alternating exact minimisation (for fixed w the ratio is a generalised Rayleigh
   quotient in w', and conversely) from many random starts at every site.  This is a numerical check, not a certificate."""
import sys
import numpy as np
from werner import d, applyB, q, q_direct, contract

t = 0.5; m = 3
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 200
rng = np.random.default_rng(2)
C = np.load(__file__.replace('code', 'data').replace('check_sci_counterexample.py', 'sci_counterexample_n3.npy'))


def gen_min(Q, H):
    ev, U = np.linalg.eigh(H)
    keep = ev > 1e-12 * ev.max()
    V = U[:, keep] / np.sqrt(ev[keep])
    lam, Wv = np.linalg.eigh(V.conj().T @ Q @ V)
    return lam[0], V @ Wv[:, 0]


def forms(Xs):
    BX = [applyB(X, m - 1, t) for X in Xs]
    Q = np.array([[np.vdot(Xs[i], BX[j]) for j in range(d)] for i in range(d)])
    H = np.array([[np.vdot(Xs[i], Xs[j]) for j in range(d)] for i in range(d)])
    return (Q + Q.conj().T) / 2, (H + H.conj().T) / 2


def min_site(s):
    best = np.inf
    for k in range(NS):
        w = rng.standard_normal(d) + 1j * rng.standard_normal(d); w /= np.linalg.norm(w)
        val = np.inf
        for it in range(300):
            lam, v = gen_min(*forms([contract(C, m, s, w, np.eye(d)[j]) for j in range(d)]))
            wp = v / np.linalg.norm(v)
            lam2, u = gen_min(*forms([contract(C, m, s, np.eye(d)[i], wp) for i in range(d)]))
            w = np.conj(u) / np.linalg.norm(u)
            if abs(lam2 - val) < 1e-15: break
            val = lam2
        X = contract(C, m, s, w, wp)
        best = min(best, q_direct(X, m - 1, t) / np.vdot(X, X).real)
    return best


nC = np.vdot(C, C).real
sv = np.linalg.svd(C, compute_uv=False)
print(f"rank check: singular values / ||C|| = {np.round(sv[:4] / np.sqrt(nC), 6)}")
R3 = q_direct(C, m, t) / nC
mins = [min_site(s) for s in range(m)]
print(f"R_3(C) = {R3:.6f}  (defining sum; via B: {q(C, m, t) / nC:.6f})")
print(f"min over contractions of R_2, per site ({NS} starts each): {np.round(mins, 6)}")
print(f"R_3(C) - (1/2) min = {R3 - 0.5 * min(mins):+.3e}   (negative: the inequality fails at this C)")
