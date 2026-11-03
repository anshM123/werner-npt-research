"""Stage 2: the energy-aware recursion of star19, checked with independent code.
(a) frame identity sum_j w_j P_u (x) P_v = 2I - F;
(b) Gamma_bar(X) = sum_j w_j Gamma_{m-1}(X_j), Sigma_m(X) = sum_j w_j Sigma_{m-1}(X_j);
(c) cost moment sum_j w_j p_j q_j = 2pq - Tr(Pc Qc) <= 2pq with Pc, Qc PSD, traces p, q;
(d) the Minkowski step: sqrt(y)(sqrt(Gbar)-sqrt(Sig_m)) <= sqrt(sum_j w_j y (sqrt G_j - sqrt S_j)_+^2)
    <= K* sqrt(sum_j w_j p_j q_j), K* = max_j atom ratio."""
import math
import torch
from core import *

g = torch.Generator().manual_seed(777)
out = []


def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)


def assemble(RA, RB, N):
    k = RA.shape[1]
    R = torch.zeros(3, N + 1, k, dtype=CD)
    R[:, :N, :] = RA.reshape(3, N, k)
    R[:, N, :] = RB
    return R.reshape(3 * (N + 1), k)


frame = octahedral_frame()
Fsw = torch.zeros(4, 4, dtype=CD)
for a in range(2):
    for b in range(2):
        Fsw[2 * a + b, 2 * b + a] = 1
S = sum(w * torch.kron(torch.outer(u, u.conj()), torch.outer(v, v.conj())) for (w, u, v) in frame)
log(f"(a) frame: ||sum w P_u(x)P_v - (2I-F)|| = {float(torch.linalg.norm(S - (2*torch.eye(4, dtype=CD) - Fsw))):.2e}, total weight {sum(w for w,_,_ in frame)}")


def qt_bil_P(P1, P2, N):
    a = P1.reshape(3, N, 3); b = P2.reshape(3, N, 3)
    return (P1.conj() * P2).sum() - 0.5 * (torch.einsum('iai->a', a).conj() * torch.einsum('iai->a', b)).sum()


def qt_bil_Q(Q1, Q2, N):
    a = Q1.reshape(3, 3, N); b = Q2.reshape(3, 3, N)
    return (Q1.conj() * Q2).sum() - 0.5 * (torch.einsum('iib->b', a).conj() * torch.einsum('iib->b', b)).sum()


def check(R, L, m, verbose=False):
    N = 2 ** m; D = N // 2; Np = N // 2
    RA, RB, LA, LB = blocks_from_factors(R, L, N)
    X, P, Q, Y = blocks(R, L, N)
    y = q1(Y); p = qtilde_P(P, N); q = qtilde_Q(Q, N)
    Gb = Gbar_m(X, D); Sm = Sigma_full(X, m)
    sG = 0.0; sS = 0.0; mom = 0.0; vec2 = 0.0; Kstar = 0.0
    for (w, u, v) in frame:
        RAu = contract_new(RA, u, D); LAv = contract_new(LA, v, D)
        Rj = assemble(RAu, RB, Np); Lj = assemble(LAv, LB, Np)
        dj = star_parts(Rj, Lj, m - 1)
        sG = sG + w * dj['Gam']; sS = sS + w * dj['Sig']
        mom = mom + w * dj['p'] * dj['q']
        a = math.sqrt(max(float(y), 0) * max(float(dj['Gam']), 0)); b = math.sqrt(max(float(y), 0) * max(float(dj['Sig']), 0))
        vec2 += w * max(a - b, 0.0) ** 2
        if float(dj['p'] * dj['q']) > 1e-28:
            Kstar = max(Kstar, float(star_ratio_from_parts(dj)))
    # cost matrices
    Pb = [contract_new(RA, torch.eye(2, dtype=CD)[a], D) @ LB.conj().T for a in range(2)]
    Qb = [RB @ contract_new(LA, torch.eye(2, dtype=CD)[a], D).conj().T for a in range(2)]
    M = torch.zeros(2, 2, dtype=CD); Nn = torch.zeros(2, 2, dtype=CD)
    for a in range(2):
        for b in range(2):
            M[b, a] = qt_bil_P(Pb[a], Pb[b], Np)
            Nn[a, b] = qt_bil_Q(Qb[a], Qb[b], Np)
    mom_formula = 2 * p * q - torch.trace(M @ Nn).real
    evM = torch.linalg.eigvalsh(0.5 * (M + M.conj().T)).min(); evN = torch.linalg.eigvalsh(0.5 * (Nn + Nn.conj().T)).min()
    lhs = math.sqrt(max(float(y), 0)) * (math.sqrt(max(float(Gb), 0)) - math.sqrt(max(float(Sm), 0)))
    mid = math.sqrt(vec2)
    rhs = Kstar * math.sqrt(max(float(mom), 0))
    sc = float(fro2(R) * fro2(L))
    return dict(eG=abs(float(sG - Gb)) / sc, eS=abs(float(sS - Sm)) / sc, eMom=abs(float(mom - mom_formula)) / sc ** 2,
                eTr=(abs(float(torch.trace(M).real - p)) + abs(float(torch.trace(Nn).real - q))) / sc,
                minev=min(float(evM), float(evN)) / sc, slack1=(mid - lhs) / sc, slack2=(rhs - mid) / sc,
                mom_le=(2 * float(p * q) - float(mom)) / sc ** 2, lhs=lhs, Kstar=Kstar)


for m in (2, 3, 4):
    N = 2 ** m
    worst = {}
    for t in range(40):
        R = rand_c(3 * (N + 1), 2, gen=g).to(CD); L = rand_c(3 * (N + 1), 2, gen=g).to(CD)
        if t % 4 == 1:   # structured: concentrate on a few strings to get positive gaps
            R[torch.rand(R.shape, generator=g) < 0.8] = 0; L[torch.rand(L.shape, generator=g) < 0.8] = 0
        r = check(R, L, m)
        for k, v in r.items():
            if k in ('slack1', 'slack2', 'mom_le', 'minev'):
                worst[k] = min(worst.get(k, 1e9), v)
            elif k not in ('lhs', 'Kstar'):
                worst[k] = max(worst.get(k, 0), v)
    log(f"m={m}: (b) |Gbar-sum|={worst['eG']:.1e} |Sig-sum|={worst['eS']:.1e}; (c) |moment-formula|={worst['eMom']:.1e} "
        f"|Tr Pc - p|+|Tr Qc - q|={worst['eTr']:.1e} min eig(Pc,Qc)={worst['minev']:.1e} min(2pq-moment)={worst['mom_le']:.1e}; "
        f"(d) min slack Minkowski={worst['slack1']:.1e} min slack K*={worst['slack2']:.1e}")

# the explicit ratio-4 operator at m=3 (positive gap example)
m = 3; N = 8
def vec(free, a):
    v = torch.zeros(3, N + 1, dtype=CD); v[free, a] = 1.0; return v.reshape(-1)
e = N; s110 = int('110', 2)
R1 = vec(0, 0) + vec(1, e); L1 = vec(1, 0) + vec(0, e); R2 = -vec(0, s110) + vec(1, e); L2 = vec(1, s110) + vec(0, e)
r = check(torch.stack([L1, L2], 1), torch.stack([R1, R2], 1), m)
log(f"explicit operator m=3 through recursion: lhs={r['lhs']:.6f}, atom K*={r['Kstar']:.6f}, slacks {r['slack1']:.2e}, {r['slack2']:.2e}")
open('stage2_recursion.log', 'w').write('\n'.join(out) + '\n')
