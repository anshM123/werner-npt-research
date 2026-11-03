"""Stage 3d: the two gauge conversions used to turn the functional bound into MB_B.
rank Y = 1 (induction17 sec 2 / orb28 sec 1): whiten R_B, rotate L_B = [v,0], free unitary -> R_B=[e0,e1];
   check Y=e0 v*, P=x1 v*, Q=e0 y1*+e1 y2*, q~(Q)=beta^2, q~(P)>=|v|^2 alpha^2/2, q1(Y)<=|v|^2,
   and the chain  sqrt(y) gap <= 16.3 |v| alpha beta <= 16.3 sqrt2 sqrt(pq).
rank Y = 2 (orb28 sec 2-3): Schur gauge Y=U diag(1,t)V*, R_B=U, L_B=V diag(1,t), R=R_A diag(1,t), L=L_A;
   check P=RV*, Q=UL*, X=RML*, p=||r||^2, q=||delta||^2, q1(Y+)>=s^2 y/2, and in the large-amplitude regime
   the denominator (7) sqrt(G)+sqrt(Gbar) >= s sqrt(y) A D d_tau and the final scaled-gap bound."""
import math
import numpy as np
import torch
from core import *

g = torch.Generator().manual_seed(99)
out = []
SQ2 = math.sqrt(2)


def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)


def gapf(X, D):
    G = float(Gam_m(X, D)); B = float(Gbar_m(X, D))
    return math.sqrt(max(G, 0)) - math.sqrt(max(B, 0)), G, B


def free_unitary(M):
    """3x3 unitary whose first two columns span... returns U3 with U3^* M = [e0 e1] for an isometry M (3x2)."""
    Qm, _ = torch.linalg.qr(torch.cat([M, rand_c(3, 1, gen=g).to(CD)], 1))
    # fix phases so that Qm[:, :2] = M exactly
    Qm = torch.cat([M, Qm[:, 2:3]], 1)
    return Qm


# ---------------- rank Y = 1
worst = dict(eY=0, eP=0, eQ=0, eq=0, p_lo=1e9, y_hi=-1e9, func=-1e9, chain=-1e9)
for D in (2, 3, 4):
    N = 2 * D
    for t in range(40):
        RA = rand_c(3 * N, 2, gen=g).to(CD); LA = rand_c(3 * N, 2, gen=g).to(CD)
        RB = rand_c(3, 2, gen=g).to(CD)
        LB = torch.outer(rand_c(3, gen=g).to(CD), rand_c(2, gen=g).to(CD).conj())          # rank one, general position
        X, P, Q, Y = RA @ LA.conj().T, RA @ LB.conj().T, RB @ LA.conj().T, RB @ LB.conj().T
        y0, p0, q0 = float(q1(Y)), float(qtilde_P(P, N)), float(qtilde_Q(Q, N))
        gap0, _, _ = gapf(X, D)
        # gauge: whiten R_B
        w, Vv = torch.linalg.eigh(RB.conj().T @ RB)
        Gm = Vv @ torch.diag((1 / torch.sqrt(w)).to(CD)) @ Vv.conj().T
        RB1, LB1, RA1, LA1 = RB @ Gm, LB @ torch.linalg.inv(Gm).conj().T, RA @ Gm, LA @ torch.linalg.inv(Gm).conj().T
        # internal unitary putting L_B = [v, 0]
        Ul, Sl, Vh = torch.linalg.svd(LB1)
        Wu = Vh.conj().T                      # L_B1 Wu = Ul diag(S) -> second column ~ 0
        RB2, LB2, RA2, LA2 = RB1 @ Wu, LB1 @ Wu, RA1 @ Wu, LA1 @ Wu
        # free unitary mapping R_B2 columns to e0, e1 (acts on the free index of every block)
        U3 = free_unitary(RB2)
        Uf = U3.conj().T
        RB3 = Uf @ RB2; LB3 = Uf @ LB2
        RA3 = kron2(Uf.contiguous(), torch.eye(N, dtype=CD)) @ RA2; LA3 = kron2(Uf.contiguous(), torch.eye(N, dtype=CD)) @ LA2
        v = LB3[:, 0]; x1, x2 = RA3[:, 0], RA3[:, 1]; y1, y2 = LA3[:, 0], LA3[:, 1]
        X3, P3, Q3, Y3 = RA3 @ LA3.conj().T, RA3 @ LB3.conj().T, RB3 @ LA3.conj().T, RB3 @ LB3.conj().T
        e0 = torch.zeros(3, dtype=CD); e0[0] = 1; e1 = torch.zeros(3, dtype=CD); e1[1] = 1
        worst['eY'] = max(worst['eY'], float(torch.linalg.norm(Y3 - torch.outer(e0, v.conj()))) + abs(float(q1(Y3)) - y0))
        worst['eP'] = max(worst['eP'], float(torch.linalg.norm(P3 - torch.outer(x1, v.conj()))) + abs(float(qtilde_P(P3, N)) - p0))
        worst['eQ'] = max(worst['eQ'], float(torch.linalg.norm(Q3 - torch.outer(e0, y1.conj()) - torch.outer(e1, y2.conj()))) + abs(float(qtilde_Q(Q3, N)) - q0))
        y1r, y2r = y1.reshape(3, N), y2.reshape(3, N)
        c = (y1r[0] + y2r[1]) / 2
        d1 = y1 - torch.kron(e0, c); d2 = y2 - torch.kron(e1, c)
        beta2 = float(fro2(d1) + fro2(d2)); alpha = math.sqrt(float(fro2(x1))); nv2 = float(fro2(v))
        worst['eq'] = max(worst['eq'], abs(q0 - beta2) / (q0 + 1e-30))
        worst['p_lo'] = min(worst['p_lo'], p0 / (nv2 * alpha ** 2 / 2))
        worst['y_hi'] = max(worst['y_hi'], y0 / nv2)
        gap3, _, _ = gapf(X3, D)
        assert abs(gap3 - gap0) < 1e-9 * (1 + abs(gap0))
        worst['func'] = max(worst['func'], gap0 / (alpha * math.sqrt(beta2)))
        worst['chain'] = max(worst['chain'], math.sqrt(y0) * gap0 / math.sqrt(p0 * q0))
log("rank Y = 1 gauge: |Y-e0v*|+|dq1|=%.1e |P-x1v*|+|dp|=%.1e |Q-..|+|dq|=%.1e |q~Q-beta^2|/q=%.1e  min p/(|v|^2a^2/2)=%.4f (>=1) "
    "max q1Y/|v|^2=%.4f (<=1)  max gap/(alpha beta)=%.4f (<=16.3)  max bridge ratio=%.4f (<=23.05)" %
    (worst['eY'], worst['eP'], worst['eQ'], worst['eq'], worst['p_lo'], worst['y_hi'], worst['func'], worst['chain']))

# ---------------- rank Y = 2
tau = 25 / 4
dtau = SQ2 - (SQ2 + 2) * (2 / tau + SQ2 / tau ** 2)
worst = dict(eP=0, eQ=0, eX=0, ep=0, eq=0, adj=1e9, den=1e9, nlarge=0, scaled=-1e9, small=-1e9, Cb=-1e9)
for D in (2, 3):
    N = 2 * D
    for t in range(60):
        if t % 2 == 0:      # generic data
            RA = rand_c(3 * N, 2, gen=g).to(CD); LA = rand_c(3 * N, 2, gen=g).to(CD)
            RB = rand_c(3, 2, gen=g).to(CD); LB = rand_c(3, 2, gen=g).to(CD)
        else:               # large scalar amplitudes, built in the canonical gauge then scrambled by a random internal G
            Vp = rand_isometry(3, 2, gen=g); Up = rand_isometry(3, 2, gen=g); t0 = 0.05 + 0.95 * float(torch.rand(1, generator=g))
            a0 = rand_c(N, gen=g).to(CD); c0 = rand_c(N, gen=g).to(CD)
            eps = 10 ** (-1 - 2 * float(torch.rand(1, generator=g)))
            R0 = torch.einsum('ik,z->izk', Vp, a0).reshape(3 * N, 2) + eps * rand_c(3 * N, 2, gen=g).to(CD)
            L0 = torch.einsum('ik,z->izk', Up, c0).reshape(3 * N, 2) + eps * rand_c(3 * N, 2, gen=g).to(CD)
            Dg0 = torch.diag(torch.tensor([1.0, t0], dtype=CD)); M0 = torch.diag(torch.tensor([1.0, 1 / t0], dtype=CD))
            Gs = rand_c(2, 2, gen=g).to(CD)
            RA, LA, RB, LB = R0 @ M0 @ Gs, L0 @ torch.linalg.inv(Gs).conj().T, Up @ Gs, Vp @ Dg0 @ torch.linalg.inv(Gs).conj().T
        X, P, Q, Y = RA @ LA.conj().T, RA @ LB.conj().T, RB @ LA.conj().T, RB @ LB.conj().T
        sc = float(torch.linalg.svdvals(Y)[0])
        RA, RB = RA / math.sqrt(sc), RB / math.sqrt(sc); LA, LB = LA / math.sqrt(sc), LB / math.sqrt(sc)
        X, P, Q, Y = RA @ LA.conj().T, RA @ LB.conj().T, RB @ LA.conj().T, RB @ LB.conj().T
        Uy, Sy, Vyh = torch.linalg.svd(Y)
        U = Uy[:, :2]; V = Vyh.conj().T[:, :2]; tt = float(Sy[1]); s = 1 / tt
        Dg = torch.diag(torch.tensor([1.0, tt], dtype=CD)); M = torch.diag(torch.tensor([1.0, s], dtype=CD))
        # gauge G with R_B G = U:  G = (U^* R_B)^{-1}
        Gm = torch.linalg.inv(U.conj().T @ RB)
        RBg, LBg, RAg, LAg = RB @ Gm, LB @ torch.linalg.inv(Gm).conj().T, RA @ Gm, LA @ torch.linalg.inv(Gm).conj().T
        assert torch.linalg.norm(RBg - U) < 1e-9 and torch.linalg.norm(LBg - V @ Dg) < 1e-8
        R = RAg @ Dg; L = LAg
        worst['eP'] = max(worst['eP'], float(torch.linalg.norm(P - R @ V.conj().T)))
        worst['eQ'] = max(worst['eQ'], float(torch.linalg.norm(Q - U @ L.conj().T)))
        worst['eX'] = max(worst['eX'], float(torch.linalg.norm(X - R @ M @ L.conj().T)))
        # scalar projections
        Rr = R.reshape(3, N, 2); a = 0.5 * torch.einsum('ik,izk->z', V.conj(), Rr)
        r = R - torch.einsum('ik,z->izk', V, a).reshape(3 * N, 2)
        Lr = L.reshape(3, N, 2); cc = 0.5 * torch.einsum('ik,izk->z', U.conj(), Lr)
        de = L - torch.einsum('ik,z->izk', U, cc).reshape(3 * N, 2)
        p, q, y = float(qtilde_P(P, N)), float(qtilde_Q(Q, N)), float(q1(Y))
        worst['ep'] = max(worst['ep'], abs(p - float(fro2(r))) / p); worst['eq'] = max(worst['eq'], abs(q - float(fro2(de))) / q)
        Yp = V @ M @ U.conj().T
        worst['adj'] = min(worst['adj'], float(q1(Yp)) / (s ** 2 * y / 2))
        A_, D_, al, be = math.sqrt(float(fro2(a))), math.sqrt(float(fro2(cc))), math.sqrt(float(fro2(r))), math.sqrt(float(fro2(de)))
        gap, G, Bb = gapf(X, D)
        if A_ * math.sqrt(y) <= tau * al or D_ * math.sqrt(y) <= tau * be:
            worst['small'] = max(worst['small'], math.sqrt(y) * gap / (al * be))
        else:
            worst['nlarge'] += 1
            worst['den'] = min(worst['den'], (math.sqrt(G) + math.sqrt(Bb)) / (s * math.sqrt(y) * A_ * D_ * dtau))
            worst['scaled'] = max(worst['scaled'], math.sqrt(y) * gap / (al * be))
        # concurrence bound and determinant bound on this datum
        Bm = internal_B_fast(R, L, D)
        Cb = concurrence(Bm)
        worst['Cb'] = max(worst['Cb'], Cb / (al * be * (21 * A_ * D_ + 4 * SQ2 * (D_ * al + A_ * be) + al * be)))
        Hval = float(H_m(X, D))
        assert Hval >= -s * Cb - 1e-9 * float(fro2(X)), (Hval, -s * Cb)
log("rank Y = 2 gauge: |P-RV*|=%.1e |Q-UL*|=%.1e |X-RML*|=%.1e |p-||r||^2|/p=%.1e |q-||d||^2|/q=%.1e  min q1(Y+)/(s^2y/2)=%.4f (>=1)  "
    "small-case max sqrt(y)gap/(ab)=%.3f (<=124.93)  large-case: n=%d min denom ratio=%s (>=1) max sqrt(y)gap/(ab)=%s  max C(B)/bound=%.4f (<=1); H>=-sC(B) all OK" %
    (worst['eP'], worst['eQ'], worst['eX'], worst['ep'], worst['eq'], worst['adj'], worst['small'], worst['nlarge'],
     ('%.3f' % worst['den']) if worst['nlarge'] else 'n/a', ('%.3f' % worst['scaled']) if worst['nlarge'] else 'n/a', worst['Cb']))
open('stage3d_gauges.log', 'w').write('\n'.join(out) + '\n')
