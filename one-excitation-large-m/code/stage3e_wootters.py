"""Stage 3e: cross-check of orb24, star25 (3) and C(B012)<=13 by maximising Wootters' concurrence DIRECTLY
(independent of the exactness of the determinant-normalised dual).  Usage: python -E stage3e_wootters.py LEMMA restarts steps seed"""
import sys, json, math, time
import torch
from core import *
from lemmas import cplx, unit, nrm, iso32, kronUc
from lemmas_conc import phys_fields, proj_tr

torch.set_num_threads(1)
SYYt = torch.tensor(SYY, dtype=CD)


def wootters_t(B):
    B = 0.5 * (B + B.conj().T)
    w, V = torch.linalg.eigh(B)
    sB = (V * torch.sqrt(torch.clamp(w, min=1e-30)).to(CD)) @ V.conj().T
    lam = torch.linalg.svdvals(sB @ SYYt @ sB.conj())
    return lam[0] - lam[1] - lam[2] - lam[3]


def make(lemma, D=2):
    if lemma == 'orb24':
        n = 6
        def f(p):
            a, p = cplx(p, (n, 2)); b, p = cplx(p, (n, 2)); c, p = cplx(p, (n, 2)); d, p = cplx(p, (n, 2))
            Vr, p = cplx(p, (n, n)); V = iso32(Vr)
            Z1 = torch.einsum('ri,rj->rij', a, b).reshape(n, 4)
            Z2 = V @ torch.einsum('ri,rj->rij', c, d).reshape(n, 4)
            Z = Z1 + Z2
            return wootters_t(Z.conj().T @ Z) / (nrm(Z1) * nrm(Z2))
        return f, 2 * (4 * n * 2 + n * n), 2.0
    if lemma == 'free13':
        Wf = torch.eye(9, dtype=CD) - 0.5 * torch.eye(9, dtype=CD)[[3 * (k % 3) + k // 3 for k in range(9)]]
        def f(p):
            s, p = p[:3], p[3:]
            Vr, p = cplx(p, (3, 2)); Ur, p = cplx(p, (3, 2)); Er, p = cplx(p, (3, 2)); Fr, p = cplx(p, (3, 2))
            ph = p[0]
            V, U = iso32(Vr), iso32(Ur)
            E = Er - 0.5 * torch.trace(V.conj().T @ Er) * V; F = Fr - 0.5 * torch.trace(U.conj().T @ Fr) * U
            E = torch.exp(s[0]) * unit(E); F = torch.exp(s[1]) * unit(F)
            al = torch.exp(s[2]) * torch.exp(1j * ph)
            Z = al * kron2(V, U) + kron2(E, U) + kron2(V, F)
            return wootters_t(Z.conj().T @ Wf @ Z) / (nrm(E) * nrm(F))
        return f, 3 + 2 * 24 + 1, 6.5
    if lemma == 'B012':
        def f(p):
            d, _ = phys_fields(p, D)
            B = gram_sum_fast([(d['R0'], d['L0']), (d['r'], d['L0']), (d['R0'], d['de'])], D)
            return wootters_t(B) / (d['A'] * d['C0'] * d['al'] * d['be'])
        return f, 4 + 2 * (6 + 6 + 2 * D + 2 * D + 12 * D + 12 * D), 13.0
    raise ValueError(lemma)


if __name__ == '__main__':
    lemma, restarts, steps, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    f, n, bound = make(lemma)
    g = torch.Generator().manual_seed(seed)
    best = -1e300; t0 = time.time()
    for r in range(restarts):
        p = torch.randn(n, generator=g)
        if lemma == 'B012':
            p[:4] = 2.0 * torch.randn(4, generator=g)
        if lemma == 'free13':
            p[:3] = 2.0 * torch.randn(3, generator=g)
        p.requires_grad_(True)
        opt = torch.optim.Adam([p], lr=0.02)
        for it in range(steps):
            val = f(p)
            if not torch.isfinite(val):
                break
            best = max(best, float(val.detach()))
            opt.zero_grad(); (-val).backward(); opt.step()
    res = dict(lemma='W_' + lemma, D=2, restarts=restarts, steps=steps, seed=seed, bound=bound, sup_found=best,
               holds=best <= bound + 1e-9, secs=time.time() - t0)
    print(json.dumps(res))
    json.dump(res, open(f'results/conc_W{lemma}_D2_s{seed}.json', 'w'))
