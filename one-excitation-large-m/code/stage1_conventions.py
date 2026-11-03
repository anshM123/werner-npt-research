"""Stage 1: conventions.  Identity of 'One free copy', Proposition M, literal vs replica forms,
explicit ratio-4 operator, exact m=3 ratio-4.018 operator (vectors copied from the published
two-sector-class/code/check_star_constant.py; evaluated with this independent code)."""
import math, sys
from fractions import Fraction
import torch
from core import *

g = torch.Generator().manual_seed(12345)
out = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    out.append(s)


def assemble(RA, RB, N):
    R = torch.zeros(3, N + 1, RA.shape[1], dtype=CD)
    R[:, :N, :] = RA.reshape(3, N, -1)
    R[:, N, :] = RB
    return R.reshape(3 * (N + 1), -1)


def B1(Z, W):
    return torch.trace(Z.conj().T @ W) - 0.5 * torch.trace(Z).conj() * torch.trace(W)


for m in (1, 2, 3):
    N = 2 ** m
    maxerr_id, maxerr_M, maxerr_lit = 0.0, 0.0, 0.0
    for trial in range(6):
        rank = 2 if trial < 3 else 5
        R = rand_c(3 * (N + 1), rank, gen=g).to(CD)
        L = rand_c(3 * (N + 1), rank, gen=g).to(CD)
        C = R @ L.conj().T
        X, P, Q, Y = blocks(R, L, N)
        lhs = (2 ** m) * qn_big(C, m)
        Xr = X.reshape(3, N, 3, N)
        X0 = torch.einsum('iaja->ij', Xr)
        rhs = Sigma_full(X, m) + q1(Y) + (2 ** m) * (qtilde_P(P, N) + qtilde_Q(Q, N)) \
            + 2 * (-1) ** m * B1(X0, Y).real
        maxerr_id = max(maxerr_id, abs(float(lhs - rhs)) / float(fro2(C)))
        # literal forms
        e1 = abs(float(Sigma_full_literal(X, m) - Sigma_full(X, m)))
        e2 = abs(float(Gamma_X_literal(X, m) - Gamma_X(X, m)))
        maxerr_lit = max(maxerr_lit, (e1 + e2) / float(fro2(C)))
        # Proposition M
        rho, tau = float(torch.rand(1, generator=g)) * 3 + 0.1, float(torch.rand(1, generator=g)) * 3 + 0.1
        d1 = torch.ones(N + 1, dtype=CD); d1[N] = (-1) ** (m + 1) * math.sqrt(rho / tau)
        d2 = torch.ones(N + 1, dtype=CD); d2[N] = math.sqrt(rho * tau)
        D1 = torch.kron(torch.eye(3, dtype=CD), torch.diag(d1))
        D2 = torch.kron(torch.eye(3, dtype=CD), torch.diag(d2))
        Z = D1 @ C @ D2
        lhsM = (2 ** m) * qn_big(C, m) - (2 / rho) * qform(Z, (3, N + 1))
        rhsM = Sigma_full(X, m) - Gamma_X(X, m) / rho + (1 - rho) * q1(Y) \
            + (2 ** m - 2 * tau) * qtilde_P(P, N) + (2 ** m - 2 / tau) * qtilde_Q(Q, N)
        maxerr_M = max(maxerr_M, abs(float(lhsM - rhsM)) / float(fro2(C)))
    log(f"m={m}: identity err {maxerr_id:.2e}; Prop M err {maxerr_M:.2e}; literal-vs-replica Sigma/Gamma err {maxerr_lit:.2e}")

# sign control: wrong sign of the cross term should fail
m, N = 2, 4
R = rand_c(3 * (N + 1), 2, gen=g).to(CD); L = rand_c(3 * (N + 1), 2, gen=g).to(CD)
C = R @ L.conj().T; X, P, Q, Y = blocks(R, L, N)
X0 = torch.einsum('iaja->ij', X.reshape(3, N, 3, N))
bad = Sigma_full(X, m) + q1(Y) + (2 ** m) * (qtilde_P(P, N) + qtilde_Q(Q, N)) - 2 * (-1) ** m * B1(X0, Y).real
log(f"control (flipped cross-term sign) err {abs(float((2**m)*qn_big(C,m)-bad))/float(fro2(C)):.2e}  (should be O(1))")

# ---- explicit ratio-4 operator of the README at m = 2..6
for m in range(2, 7):
    N = 2 ** m
    def vec(free, a):
        v = torch.zeros(3, N + 1, dtype=CD); v[free, a] = 1.0; return v.reshape(-1)
    e = N
    s110 = int('11' + '0' * (m - 2), 2)
    R1 = vec(0, 0) + vec(1, e); L1 = vec(1, 0) + vec(0, e)
    R2 = -vec(0, s110) + vec(1, e); L2 = vec(1, s110) + vec(0, e)
    # C^* = |R1><L1| + |R2><L2|  =>  C = |L1><R1| + |L2><R2|
    Rf = torch.stack([L1, L2], 1); Lf = torch.stack([R1, R2], 1)
    d = star_parts(Rf, Lf, m)
    log(f"explicit operator m={m}: Sigma={float(d['Sig']):.6f} Gamma={float(d['Gam']):.6f} q1Y={float(d['y']):.6f} "
        f"p={float(d['p']):.6f} q={float(d['q']):.6f} ratio={float(star_ratio_from_parts(d)):.10f}")

# ---- exact m=3 operator with ratio 4.018 (README "4 is not the supremum")
A_cols = [(324+868j, -1+1j), (177-115j, 2-3j), (345+31j, -4-2j), (18-54j, 4-2j), (119-26j, -2j), (-44-4j, 5j),
          (48-50j, 1-5j), (10-10j, -1j)]
B_cols = [(-807, 41), (-175+317j, -161-12j), (133+279j, -154+52j), (97+35j, -38-439j), (-66+353j, -12+50j),
          (104+63j, 574+187j), (121+12j, -421-270j), (18-39j, -406-136j)]
m, N = 3, 8
a1 = torch.tensor([complex(r[0]) for r in A_cols], dtype=CD) / 1024
a2 = torch.tensor([complex(r[1]) for r in A_cols], dtype=CD) / 1024
b1 = torch.tensor([complex(r[0]) for r in B_cols], dtype=CD) / 1024
b2 = torch.tensor([complex(r[1]) for r in B_cols], dtype=CD) / 1024


def vecA(free, a):
    v = torch.zeros(3, N + 1, dtype=CD); v[free, :N] = a; return v.reshape(-1)


def vecE(free):
    v = torch.zeros(3, N + 1, dtype=CD); v[free, N] = 1.0; return v.reshape(-1)


r1 = vecA(0, a1); l1 = vecA(1, b1); r2 = vecA(0, a2) + vecE(1); l2 = vecA(1, b2) + vecE(0)
Rf = torch.stack([r1, r2], 1); Lf = torch.stack([l1, l2], 1)
d = star_parts(Rf, Lf, m)
log(f"exact m=3 operator: Gamma*2^40={float(d['Gam'])*2**40:.1f} (README 1358550811983), "
    f"Sigma*2^40={float(d['Sig'])*2**40:.1f} (README 1308487705385), q1Y={float(d['y']):.6f}, ratio={float(star_ratio_from_parts(d)):.10f}")
torch.save({'R': Rf, 'L': Lf}, 'exact_m3_operator.pt')

open('stage1_conventions.log', 'w').write('\n'.join(out) + '\n')
