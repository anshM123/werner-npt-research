"""Stage 3b: concurrence-side lemmas of the rank-2-Y branch, plus the adjugate comparison.
Objectives use the witness value  w(M,B) = -Tr(W_M B)/|det M|, which is <= C(B) for EVERY M
(proved in the report: <psi|W_M|psi> >= -2|det M||det Psi|, then convex roof), so a value above a
claimed bound is a genuine counterexample.  Best points are re-evaluated with Wootters' formula.
Usage: python -E lemmas_conc.py LEMMA D restarts steps seed"""
import sys, json, math, time
import numpy as np
import torch
from core import *
from lemmas import cplx, unit, nrm, iso32, kronUc

torch.set_num_threads(1)
SQ2 = math.sqrt(2.0)


def wit(Mp, B):
    M = torch.complex(Mp[:4].reshape(2, 2), Mp[4:8].reshape(2, 2))
    W = torch.einsum('kl,ij->kjil', M, M.conj()).reshape(4, 4)
    return -(torch.trace(W @ B)).real / (torch.linalg.det(M).abs() + 1e-300)


def proj_tr(Rz, V):
    """make every sector slice trace-orthogonal to the plane V: R has rows (free, z); slice R_z is 3x2."""
    n = Rz.shape[0] // 3
    Rr = Rz.reshape(3, n, 2)
    t = torch.einsum('ik,izk->z', V.conj(), Rr)          # Tr(V^* R_z)
    return (Rr - 0.5 * torch.einsum('ik,z->izk', V, t)).reshape(3 * n, 2)


def phys_fields(p, D):
    """R = V (x) a + r, L = U (x) c + delta, trace-orthogonal residuals, explicit log scales."""
    s, p = p[:4], p[4:]
    Vr, p = cplx(p, (3, 2)); Ur, p = cplx(p, (3, 2))
    ar, p = cplx(p, (2 * D,)); cr, p = cplx(p, (2 * D,))
    rr, p = cplx(p, (6 * D, 2)); dr, p = cplx(p, (6 * D, 2))
    V, U = iso32(Vr), iso32(Ur)
    a = torch.exp(s[0]) * unit(ar); c = torch.exp(s[1]) * unit(cr)
    r = torch.exp(s[2]) * unit(proj_tr(rr, V)); de = torch.exp(s[3]) * unit(proj_tr(dr, U))
    return dict(V=V, U=U, a=a, c=c, r=r, de=de, A=nrm(a), C0=nrm(c), al=nrm(r), be=nrm(de),
                R0=kronUc(V, a), L0=kronUc(U, c)), p


def nphys(D):
    return 4 + 2 * (6 + 6 + 2 * D + 2 * D + 12 * D + 12 * D)


def make(lemma, D):
    if lemma == 'dual':        # dual exactness: sup_M w(M,B) vs Wootters C(B) for a fixed random PSD B
        return None
    if lemma == 'orb24':       # C(G(Z1+Z2)) <= 2||Z1|| ||Z2||, G(Z1), G(Z2) separable; env dim n
        n = 6
        def f(p):
            Mp, p = p[:8], p[8:]
            a, p = cplx(p, (n, 2)); b, p = cplx(p, (n, 2)); c, p = cplx(p, (n, 2)); d, p = cplx(p, (n, 2))
            Vr, p = cplx(p, (n, n)); V = iso32(Vr)
            Z1 = torch.einsum('ri,rj->rij', a, b).reshape(n, 4)
            Z2 = V @ torch.einsum('ri,rj->rij', c, d).reshape(n, 4)
            Z = Z1 + Z2
            return wit(Mp, Z.conj().T @ Z) / (nrm(Z1) * nrm(Z2))
        return f, 8 + 2 * (4 * n * 2 + n * n), 2.0
    if lemma == 'free13':      # star25 (3): C(Z^* W Z) <= (13/2) e f, Z = alpha V(x)U + E(x)U + V(x)F
        Wf = torch.eye(9, dtype=CD) - 0.5 * torch.eye(9, dtype=CD)[[3 * (k % 3) + k // 3 for k in range(9)]]
        def f(p):
            Mp, p = p[:8], p[8:]
            s, p = p[:3], p[3:]
            Vr, p = cplx(p, (3, 2)); Ur, p = cplx(p, (3, 2)); Er, p = cplx(p, (3, 2)); Fr, p = cplx(p, (3, 2))
            ph = p[0]
            V, U = iso32(Vr), iso32(Ur)
            E = Er - 0.5 * torch.trace(V.conj().T @ Er) * V; F = Fr - 0.5 * torch.trace(U.conj().T @ Fr) * U
            E = torch.exp(s[0]) * unit(E); F = torch.exp(s[1]) * unit(F)
            al = torch.exp(s[2]) * torch.exp(1j * ph)
            Z = al * kron2(V, U) + kron2(E, U) + kron2(V, F)
            return wit(Mp, Z.conj().T @ Wf @ Z) / (nrm(E) * nrm(F))
        return f, 8 + 3 + 2 * 24 + 1, 6.5
    if lemma == 'B012':        # induction24 (7) / star25 (8): C(B012) <= 13 A C0 alpha beta
        def f(p):
            Mp, p = p[:8], p[8:]
            d, _ = phys_fields(p, D)
            B = gram_sum_fast([(d['R0'], d['L0']), (d['r'], d['L0']), (d['R0'], d['de'])], D)
            return wit(Mp, B) / (d['A'] * d['C0'] * d['al'] * d['be'])
        return f, 8 + nphys(D), 13.0
    if lemma == 'Bfull':       # orb26 (11) as used in orb28 (5): C(B) <= ab[21AD + 4sqrt2(D a + A b) + a b]
        def f(p):
            Mp, p = p[:8], p[8:]
            d, _ = phys_fields(p, D)
            B = internal_B_fast(d['R0'] + d['r'], d['L0'] + d['de'], D)
            al, be, A, C0 = d['al'], d['be'], d['A'], d['C0']
            bound = al * be * (21 * A * C0 + 4 * SQ2 * (C0 * al + A * be) + al * be)
            return wit(Mp, B) / bound
        return f, 8 + nphys(D), 1.0
    if lemma == 'B3':          # orb26 (2): C(B3) <= ||E||^2 ||F||^2, B3 = (E^*(x)F^*)K(E(x)F)
        def f(p):
            Mp, p = p[:8], p[8:]
            E, p = cplx(p, (6 * D, 2)); F, p = cplx(p, (6 * D, 2))
            return wit(Mp, internal_B_fast(E, F, D)) / (fro2(E) * fro2(F))
        return f, 8 + 2 * 24 * D, 1.0
    if lemma == 'comp':        # star21 (7): z0<=2sqrt2 A C0, z1<=2 C0 al, z2<=2 A be, z3<=sqrt2 al be  (max of the 4 ratios)
        def f(p):
            d, _ = phys_fields(p, D)
            z = lambda R, L: torch.sqrt(torch.trace(internal_B_fast(R, L, D)).real)
            r0 = z(d['R0'], d['L0']) / (2 * SQ2 * d['A'] * d['C0'])
            r1 = z(d['r'], d['L0']) / (2 * d['C0'] * d['al'])
            r2 = z(d['R0'], d['de']) / (2 * d['A'] * d['be'])
            r3 = z(d['r'], d['de']) / (SQ2 * d['al'] * d['be'])
            return torch.max(torch.stack([r0, r1, r2, r3]))
        return f, nphys(D), 1.0
    if lemma == 'detbound':    # H(R M L^*) >= -|det M| C(B): here check H(RML^*) = Tr(W_M B) exactly (value = error)
        return None
    if lemma == 'adj':         # q1(Y^+) >= s^2 y / 2  for Y = U diag(1,t) V^*  -> maximise 1 - 2 t^2 q1(Y^+)/q1(Y)
        def f(p):
            Ur, p = cplx(p, (3, 2)); Vr, p = cplx(p, (3, 2)); t = torch.sigmoid(p[0])
            U, V = iso32(Ur), iso32(Vr)
            Y = U @ torch.diag(torch.stack([torch.tensor(1.0, dtype=CD), t.to(CD)])) @ V.conj().T
            Yp = V @ torch.diag(torch.stack([torch.tensor(1.0, dtype=CD), (1 / t).to(CD)])) @ U.conj().T
            return 1 - 2 * t ** 2 * q1(Yp) / q1(Y)
        return f, 24 + 1, 0.0
    raise ValueError(lemma)


def run(lemma, D, restarts, steps, seed):
    f, n, bound = make(lemma, D)
    g = torch.Generator().manual_seed(seed)
    best, bestp = -1e300, None
    t0 = time.time()
    for r_ in range(restarts):
        p = torch.randn(n, generator=g)
        if lemma in ('B012', 'Bfull', 'comp'):
            off = 8 if lemma != 'comp' else 0
            p[off:off + 4] = 2.5 * torch.randn(4, generator=g)
        if lemma == 'free13':
            p[8:11] = 2.0 * torch.randn(3, generator=g)
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
                if lemma in ('B012', 'Bfull', 'comp'):
                    off = 8 if lemma != 'comp' else 0
                    p[off:off + 4].clamp_(-10, 10)
                if lemma == 'free13':
                    p[8:11].clamp_(-10, 10)
    return dict(lemma=lemma, D=D, restarts=restarts, steps=steps, seed=seed, bound=bound, sup_found=best,
                holds=best <= bound + 1e-9, secs=time.time() - t0), bestp


if __name__ == '__main__':
    lemma, D, restarts, steps, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    res, bp = run(lemma, D, restarts, steps, seed)
    print(json.dumps(res))
    with open(f'results/conc_{lemma}_D{D}_s{seed}.json', 'w') as fh:
        json.dump(res, fh)
    torch.save(bp, f'results/conc_{lemma}_D{D}_s{seed}.pt')
