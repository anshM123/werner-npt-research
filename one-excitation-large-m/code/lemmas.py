"""Stage 3: adversarial tests of every lemma used by the bridge MB_125 (orb37 + orb28 + inputs).
Usage: python -E lemmas.py LEMMA D restarts steps seed
Each objective returns a scalar to MAXIMISE; the lemma holds iff the sup is <= the stated bound."""
import sys, json, math, time
import torch
from core import *

torch.set_num_threads(1)
SQ2 = math.sqrt(2.0)


def nrm(v):
    return torch.sqrt(fro2(v))


def cplx(p, shape):
    n = int(torch.tensor(shape).prod())
    return torch.complex(p[:n].reshape(shape), p[n:2 * n].reshape(shape)), p[2 * n:]


def unit(v):
    return v / (nrm(v) + 1e-300)


def svals(X):
    s = torch.linalg.svdvals(X)
    return s[0], (s[1] if s.numel() > 1 else torch.tensor(0.0))


def iso32(M):
    Qm, Rm = torch.linalg.qr(M)
    return Qm


def kronUc(U, c):
    return torch.einsum('ik,z->izk', U, c).reshape(3 * c.numel(), 2)


def FL_data(p, D):
    """F=[f1,f2], U, c, delta with explicit log scales.  Returns X=F L^*, alpha, T, S, beta and pieces."""
    s, p = p[:4], p[4:]
    f1r, p = cplx(p, (6 * D,)); f2r, p = cplx(p, (6 * D,))
    Ur, p = cplx(p, (3, 2)); cr, p = cplx(p, (2 * D,)); dr, p = cplx(p, (6 * D, 2))
    f1 = torch.exp(s[0]) * unit(f1r); f2 = torch.exp(s[1]) * unit(f2r)
    c = torch.exp(s[2]) * unit(cr); delta = torch.exp(s[3]) * unit(dr)
    U = iso32(Ur)
    F = torch.stack([f1, f2], 1)
    L0 = kronUc(U, c)
    L = L0 + delta
    X = F @ L.conj().T
    return dict(X=X, F=F, L=L, L0=L0, U=U, c=c, delta=delta, alpha=nrm(f1), T=nrm(f2), S=nrm(c), beta=nrm(delta))


def nFL(D):
    return 4 + 2 * (6 * D * 2 + 6 + 2 * D + 12 * D)


def rank2(p, d1, d2):
    a, p = cplx(p, (d1 * d2, 2)); b, p = cplx(p, (d1 * d2, 2))
    return a @ b.conj().T


def gap_of(X, D):
    G = torch.clamp(Gam_m(X, D), min=0.0); B = torch.clamp(Gbar_m(X, D), min=0.0)
    return (G - B) / (torch.sqrt(G) + torch.sqrt(B) + 1e-300), G, B


# ---------------------------------------------------------------- objectives
def make(lemma, D):
    if lemma == 'flat2':      # induction16 (4): q2(Z) >= (s1-s2)^2/4 on d1 x d2 (here 3 x 2D);  maximise deficit
        n = 2 * (2 * 6 * D * 2)
        def f(p):
            Z = rank2(p, 3, 2 * D); s1, s2 = svals(Z)
            return ((s1 - s2) ** 2 / 4 - qform(Z, (3, 2 * D))) / fro2(Z)
        return f, n, 0.0
    if lemma == 'flat2b':     # same, two arbitrary factors 2D x 2D
        n = 2 * (2 * 4 * D * D * 2)
        def f(p):
            Z = rank2(p, 2 * D, 2 * D); s1, s2 = svals(Z)
            return ((s1 - s2) ** 2 / 4 - qform(Z, (2 * D, 2 * D))) / fro2(Z)
        return f, n, 0.0
    if lemma == 'flat3':      # orb33 (1): Gbar >= (s1-s2)^2/2
        n = 2 * (2 * 6 * D * 2)
        def f(p):
            Z = rank2(p, 3, 2 * D); s1, s2 = svals(Z)
            return ((s1 - s2) ** 2 / 2 - Gbar_m(Z, D)) / fro2(Z)
        return f, n, 0.0
    if lemma == 'Hpure_hi':   # star20 (1): H(uv^*) <= 2|u|^2|v|^2
        n = 2 * (2 * 6 * D)
        def f(p):
            u, p = cplx(p, (6 * D,)); v, p = cplx(p, (6 * D,))
            return H_m(torch.outer(u, v.conj()), D) / (fro2(u) * fro2(v))
        return f, n, 2.0
    if lemma == 'Hpure_lo':   # H(uv^*) >= 0
        n = 2 * (2 * 6 * D)
        def f(p):
            u, p = cplx(p, (6 * D,)); v, p = cplx(p, (6 * D,))
            return -H_m(torch.outer(u, v.conj()), D) / (fro2(u) * fro2(v))
        return f, n, 0.0
    if lemma == 'Hcurv':      # star20 (6): H >= -4 s1 s2
        n = 2 * (2 * 6 * D * 2)
        def f(p):
            Z = rank2(p, 3, 2 * D); s1, s2 = svals(Z)
            return -(H_m(Z, D) + 4 * s1 * s2) / fro2(Z)
        return f, n, 0.0
    if lemma == 'Gup':        # Gamma <= 2||Z||^2 (rank 2)
        n = 2 * (2 * 6 * D * 2)
        def f(p):
            Z = rank2(p, 3, 2 * D); return Gam_m(Z, D) / fro2(Z)
        return f, n, 2.0
    if lemma == 'Gbup':       # Gamma_bar <= 4||Z||^2
        n = 2 * (2 * 6 * D * 2)
        def f(p):
            Z = rank2(p, 3, 2 * D); return Gbar_m(Z, D) / fro2(Z)
        return f, n, 4.0
    if lemma == 'Gpos':       # Gamma >= 0 and Gamma_bar >= 0 (rank 2): maximise -min
        n = 2 * (2 * 6 * D * 2)
        def f(p):
            Z = rank2(p, 3, 2 * D); return -torch.minimum(Gam_m(Z, D), Gbar_m(Z, D)) / fro2(Z)
        return f, n, 0.0
    if lemma == 'safe':       # induction17 (5): H(Z) >= 0 for co-range in U (x) c
        n = 2 * (6 * D * 2 + 6 + 2 * D)
        def f(p):
            F, p = cplx(p, (6 * D, 2)); Ur, p = cplx(p, (3, 2)); c, p = cplx(p, (2 * D,))
            Z = F @ kronUc(iso32(Ur), c).conj().T
            return -H_m(Z, D) / fro2(Z)
        return f, n, 0.0
    if lemma == 'aux':        # orb33 (10): sqrt G - sqrt Gbar <= 2 sqrt2 s2
        n = 2 * (2 * 6 * D * 2)
        def f(p):
            Z = rank2(p, 3, 2 * D); s1, s2 = svals(Z); gp, _, _ = gap_of(Z, D)
            return gp / (s2 + 1e-300)
        return f, n, 2 * SQ2
    if lemma == 'curv27':     # orb27 (4): -H(FL^*) <= 8 sqrt2 a b T S + 4 a T b^2
        def f(p):
            d = FL_data(p, D)
            bound = 8 * SQ2 * d['alpha'] * d['beta'] * d['T'] * d['S'] + 4 * d['alpha'] * d['T'] * d['beta'] ** 2
            return -H_m(d['X'], D) / bound
        return f, nFL(D), 1.0
    if lemma == 'func':       # orb37 (10): gap <= 163/10 alpha beta
        def f(p):
            d = FL_data(p, D); gp, _, _ = gap_of(d['X'], D)
            return gp / (d['alpha'] * d['beta'])
        return f, nFL(D), 16.3
    if lemma == 'safecomp':   # orb27 sec 2 / orb37 (7): gap <= (sqrt2+2) sqrt(T^2+alpha^2) beta
        def f(p):
            d = FL_data(p, D); gp, _, _ = gap_of(d['X'], D)
            return gp / (torch.sqrt(d['T'] ** 2 + d['alpha'] ** 2) * d['beta'])
        return f, nFL(D), SQ2 + 2
    if lemma == 'denom':      # orb37 (5): sqrtG + sqrtGbar >= sqrt2 [TS - aS - b(T+a)]
        def f(p):
            d = FL_data(p, D); _, G, B = gap_of(d['X'], D)
            rhs = SQ2 * (d['T'] * d['S'] - d['alpha'] * d['S'] - d['beta'] * (d['T'] + d['alpha']))
            return (rhs - torch.sqrt(G) - torch.sqrt(B)) / (d['T'] * d['S'])
        return f, nFL(D), 0.0
    raise ValueError(lemma)


def run(lemma, D, restarts, steps, seed):
    f, n, bound = make(lemma, D)
    g = torch.Generator().manual_seed(seed)
    best, bestp = -1e300, None
    t0 = time.time()
    for r in range(restarts):
        p = torch.randn(n, generator=g)
        if lemma in ('curv27', 'func', 'safecomp', 'denom'):
            p[:4] = 3.0 * torch.randn(4, generator=g)          # explore scale regimes
        p.requires_grad_(True)
        opt = torch.optim.Adam([p], lr=0.03)
        for it in range(steps):
            opt.zero_grad()
            val = f(p)
            if not torch.isfinite(val):
                break
            v = float(val.detach())
            if v > best:
                best, bestp = v, p.detach().clone()
            (-val).backward()
            opt.step()
            with torch.no_grad():
                if lemma in ('curv27', 'func', 'safecomp', 'denom'):
                    p[:4].clamp_(-12, 12)
    return dict(lemma=lemma, D=D, restarts=restarts, steps=steps, seed=seed, bound=bound, sup_found=best,
                holds=best <= bound + 1e-9, secs=time.time() - t0), bestp


if __name__ == '__main__':
    lemma, D, restarts, steps, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    res, bp = run(lemma, D, restarts, steps, seed)
    print(json.dumps(res))
    with open(f'results/lemma_{lemma}_D{D}_s{seed}.json', 'w') as fh:
        json.dump(res, fh)
    torch.save(bp, f'results/lemma_{lemma}_D{D}_s{seed}.pt')
