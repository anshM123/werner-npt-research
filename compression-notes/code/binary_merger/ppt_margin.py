"""PPT margin of the one-atom form.
For a target X (Tr X = 1) and an atom G^ = G(E, F)/Tr G(E, F):
    mu(E, F, lam) = min( lambda_min(X - lam G^), lambda_min((X - lam G^)^Gamma) ),
    mu*(X) = max over (E, F, lam >= 0) of mu.
mu*(X) >= 0 is necessary for SF at X (the remainder must be PSD and PPT), and mu*(X) > 0 means a PPT remainder with a
margin.  The maximisation uses a soft-min of the 18 eigenvalues (beta) and L-BFGS from several starts; the reported
value is the exact (hard-min) mu at the final parameters, so it is a LOWER bound on mu*(X).
"""
import numpy as np
import torch
from bm_core import RHO_T, D, N, pt, G

torch.set_default_dtype(torch.float64)


def pt_t(Y):
    return Y.reshape(D, D, D, D).permute(0, 3, 2, 1).reshape(N, N)


def atom_hat_t(E, Fm):
    K = torch.kron(E, Fm); Y = K @ RHO_T @ K.conj().T
    Y = (Y + Y.conj().T) / 2
    return Y / torch.trace(Y).real


def softmin(ev, beta):
    return -torch.logsumexp(-beta * ev, 0) / beta


def mu_hard(X, E, Fm, lam):
    Gh = G(E, Fm); Gh = Gh / np.trace(Gh).real
    S = X - lam * Gh
    return min(np.linalg.eigvalsh(S)[0], np.linalg.eigvalsh(pt(S))[0])


def maximize_mu(X, starts, beta=3000.0, iters=300, lam0=None):
    """returns (mu_hard, E, F, lam) for the best start; X Hermitian with Tr X = 1"""
    Xt = torch.tensor(X, dtype=torch.complex128)
    best = (-np.inf, None, None, None)
    for (E0, F0) in starts:
        E = torch.tensor(E0, dtype=torch.complex128).clone().requires_grad_()
        Fm = torch.tensor(F0, dtype=torch.complex128).clone().requires_grad_()
        # initial lam: best hard margin on a grid
        Gh0 = G(E0, F0); Gh0 /= np.trace(Gh0).real
        lams = np.linspace(0, 1.5, 61)
        vals = [min(np.linalg.eigvalsh(X - l * Gh0)[0], np.linalg.eigvalsh(pt(X - l * Gh0))[0]) for l in lams]
        l0 = lams[int(np.argmax(vals))] if lam0 is None else lam0
        l0 = min(max(l0, 1e-3), 1.9); sl = torch.tensor(np.log(l0 / (2 - l0)), requires_grad=True)   # lam = 2 sigmoid(sl)
        opt = torch.optim.LBFGS([E, Fm, sl], lr=1, max_iter=iters, history_size=30, line_search_fn='strong_wolfe',
                                tolerance_grad=1e-12, tolerance_change=1e-14)

        def cl():
            opt.zero_grad()
            S = Xt - 2 * torch.sigmoid(sl) * atom_hat_t(E, Fm)
            ev = torch.cat([torch.linalg.eigvalsh(S), torch.linalg.eigvalsh(pt_t(S))])
            loss = -softmin(ev, beta)
            loss.backward(); return loss
        try:
            opt.step(cl)
        except Exception:
            pass
        En, Fn, ln = E.detach().numpy(), Fm.detach().numpy(), float(2 / (1 + np.exp(-sl.detach().numpy())))
        try:
            m = mu_hard(X, En, Fn, ln) if np.all(np.isfinite(En)) and np.all(np.isfinite(Fn)) else -np.inf
        except np.linalg.LinAlgError:
            m = -np.inf
        if m > best[0]:
            best = (m, En, Fn, ln)
    return best
