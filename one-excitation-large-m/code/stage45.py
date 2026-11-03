"""Stage 4: direct adversarial search for the merged bridge ratio
      sqrt(q1 Y)(sqrt Gamma - sqrt Gamma_bar)/sqrt(q~P q~Q)     (MB; claimed <= 125, <= 23.05 for rank Y <= 1)
   on C^3 (x) ((C^2 (x) C^D) (+) C e).
Stage 5: direct adversarial search for the (star) ratio
      sqrt(q1 Y)(sqrt Gamma_X - sqrt Sigma_full)/sqrt(q~P q~Q)  (needed <= 2^m - 2)
   at level m.
Usage: python -E stage45.py {bridge|bridge1|star|star1} PARAM restarts steps seed [init]
   PARAM = D for bridge, m for star.  init in {rand, w4, x4018} (seeding near known extremal operators)."""
import sys, json, math, time
import torch
from core import *
from lemmas import cplx, unit, nrm

torch.set_num_threads(1)


def fields(p, N, rank1Y):
    s, p = p[:4], p[4:]
    RA, p = cplx(p, (3 * N, 2)); RB, p = cplx(p, (3, 2)); LA, p = cplx(p, (3 * N, 2)); LB, p = cplx(p, (3, 2))
    if rank1Y:
        LB = torch.stack([LB[:, 0], torch.zeros(3, dtype=CD)], 1)
    RA = torch.exp(s[0]) * RA; RB = torch.exp(s[1]) * RB; LA = torch.exp(s[2]) * LA; LB = torch.exp(s[3]) * LB
    return RA, RB, LA, LB


def npar(N):
    return 4 + 2 * (2 * (3 * N * 2) + 2 * 6)


def ratio_terms(RA, RB, LA, LB, N, kind, D_or_m):
    X = RA @ LA.conj().T; P = RA @ LB.conj().T; Q = RB @ LA.conj().T; Y = RB @ LB.conj().T
    y = q1(Y); pp = qtilde_P(P, N); qq = qtilde_Q(Q, N)
    if kind.startswith('bridge'):
        G1 = Gam_m(X, D_or_m); G2 = Gbar_m(X, D_or_m)
    else:
        G1 = Gamma_X(X, D_or_m); G2 = Sigma_full(X, D_or_m)
    G1 = torch.clamp(G1, min=0.0); G2 = torch.clamp(G2, min=0.0)
    num = torch.sqrt(torch.clamp(y, min=0.0)) * (G1 - G2) / (torch.sqrt(G1) + torch.sqrt(G2) + 1e-300)
    return num / torch.sqrt(torch.clamp(pp * qq, min=1e-300)), dict(y=y, p=pp, q=qq, G1=G1, G2=G2)


def to_params(RA, RB, LA, LB):
    flat = [RA, RB, LA, LB]
    re = torch.cat([t.real.reshape(-1) for t in flat]); im = torch.cat([t.imag.reshape(-1) for t in flat])
    # cplx() reads each block as (real block, imag block) sequentially
    out = []
    for t in flat:
        out += [t.real.reshape(-1), t.imag.reshape(-1)]
    return torch.cat([torch.zeros(4)] + out)


def seed_fields(init, N, g, eps=0.05):
    """w4: the explicit ratio-4 witness (strings 0..0 and 110..0); x4018: the exact m=3 operator embedded."""
    def zeros():
        return torch.zeros(3, N, 2, dtype=CD), torch.zeros(3, 2, dtype=CD), torch.zeros(3, N, 2, dtype=CD), torch.zeros(3, 2, dtype=CD)
    RA, RB, LA, LB = zeros()
    if init == 'w4':
        m = int(round(math.log2(N)))
        s110 = int('11' + '0' * (m - 2), 2)
        # C = |L1><R1| + |L2><R2|:  range factors (L1, L2), co-range factors (R1, R2) of the README C^*
        # L1 = |1>|0..0> + |0>|e>, L2 = |1>|110..0> + |0>|e>, R1 = |0>|0..0> + |1>|e>, R2 = -|0>|110..0> + |1>|e>
        RA[1, 0, 0] = 1; RB[0, 0] = 1; RA[1, s110, 1] = 1; RB[0, 1] = 1
        LA[0, 0, 0] = 1; LB[1, 0] = 1; LA[0, s110, 1] = -1; LB[1, 1] = 1
    elif init == 'x4018':
        d = torch.load('exact_m3_operator.pt')
        R3, L3 = d['R'].reshape(3, 9, 2), d['L'].reshape(3, 9, 2)
        k = N // 8          # append product factors |0> on the extra copies: A index a -> a * k
        for a in range(8):
            RA[:, a * k, :] = R3[:, a, :]; LA[:, a * k, :] = L3[:, a, :]
        RB[:] = R3[:, 8, :]; LB[:] = L3[:, 8, :]
    RA = RA.reshape(3 * N, 2); LA = LA.reshape(3 * N, 2)
    p = to_params(RA, RB, LA, LB)
    p[4:] += eps * torch.randn(p.numel() - 4, generator=g)
    return p


def run(kind, param, restarts, steps, seed, init):
    rank1 = kind.endswith('1')
    if kind.startswith('bridge'):
        D = param; N = 2 * D; arg = D
    else:
        m = param; N = 2 ** m; arg = m
    n = npar(N)
    g = torch.Generator().manual_seed(seed)
    best, bestinfo, bestp = -1e300, None, None
    t0 = time.time()
    for r_ in range(restarts):
        if init == 'rand':
            p = torch.randn(n, generator=g); p[:4] = 1.5 * torch.randn(4, generator=g)
        else:
            p = seed_fields(init, N, g, eps=(1e-5 if init == 'x4018' else 0.02) * (1 + r_))
        p.requires_grad_(True)
        opt = torch.optim.Adam([p], lr=(0.02 if init == 'rand' else 0.002))
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps, eta_min=1e-4)
        for it in range(steps):
            opt.zero_grad()
            val, info = ratio_terms(*fields(p, N, rank1), N, kind, arg)
            if not torch.isfinite(val):
                break
            v = float(val.detach())
            if v > best:
                best = v; bestp = p.detach().clone()
                bestinfo = {k: float(t.detach()) for k, t in info.items()}
            (-val).backward()
            opt.step(); sched.step()
            with torch.no_grad():
                p[:4].clamp_(-8, 8)
    # re-evaluate best point and its rank
    with torch.no_grad():
        RA, RB, LA, LB = fields(bestp, N, rank1)
        val, info = ratio_terms(RA, RB, LA, LB, N, kind, arg)
        Yr = torch.linalg.svdvals(RB @ LB.conj().T)
    res = dict(kind=kind, param=param, init=init, restarts=restarts, steps=steps, seed=seed, sup_found=float(val),
               info={k: float(t) for k, t in info.items()}, Y_svals=[float(x) for x in Yr], secs=time.time() - t0)
    if kind.startswith('bridge'):
        res['claimed_bound'] = 23.05 if rank1 else 125.0
    else:
        res['needed'] = 2 ** param - 2
    return res, bestp


if __name__ == '__main__':
    kind, param, restarts, steps, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    init = sys.argv[6] if len(sys.argv) > 6 else 'rand'
    res, bp = run(kind, param, restarts, steps, seed, init)
    print(json.dumps(res))
    tag = f'{kind}_{param}_{init}_s{seed}'
    with open(f'results/s45_{tag}.json', 'w') as fh:
        json.dump(res, fh)
    torch.save(bp, f'results/s45_{tag}.pt')
