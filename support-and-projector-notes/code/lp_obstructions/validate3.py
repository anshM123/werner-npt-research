"""Validate the three-replica bookkeeping on actual rank-two projectors (n = 1, 2, 3)."""
import itertools, string, sys
import numpy as np
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'global_reductions'))
from qtools import q_dense, weyl, rand_c
from rep3 import (PERMS, STD, sign, ncycles, compose, inverse, DIMV, dimM, pi_j, const, sym_power_poly,
                  sign_projector, group_action)

LET = string.ascii_letters


def moments(Pi, n):
    T = Pi.reshape((3,) * (2 * n))
    out = {}
    for sig in itertools.product(PERMS, repeat=n):
        # letters r[k][nu] for row index of replica k at copy nu
        r = [[LET[3 * nu + k] for nu in range(n)] for k in range(3)]
        specs = []
        for i in range(3):
            rows = ''.join(r[i][nu] for nu in range(n))
            cols = ''.join(r[inverse(sig[nu])[i]][nu] for nu in range(n))
            specs.append(rows + cols)
        out[sig] = np.einsum(','.join(specs) + '->', T, T, T, optimize=True)
    return out


def coefficients(mom, n):
    G = np.array([[3.0 ** ncycles(compose(p, q)) for q in PERMS] for p in PERMS])
    Gi = np.linalg.inv(G)
    idx = {p: i for i, p in enumerate(PERMS)}
    c = {}
    for tau in itertools.product(PERMS, repeat=n):
        val = 0
        for sig, mv in mom.items():
            w = 1.0
            for nu in range(n):
                w *= Gi[idx[tau[nu]], idx[sig[nu]]]
            val += w * mv
        c[tau] = val
    return c


def rho(lam, p):
    if lam == 't':
        return np.array([[1.0]])
    if lam == 'a':
        return np.array([[float(sign(p))]])
    return STD[p]


def block(c, lamv):
    dim = 2 ** sum(1 for l in lamv if l == 's')
    B = np.zeros((dim, dim), complex)
    for tau, cv in c.items():
        M = np.array([[1.0]])
        for nu, l in enumerate(lamv):
            M = np.kron(M, rho(l, tau[nu]))
        B += cv * M
    return B


rng = np.random.default_rng(2)
ok = True
# character identity for pi_j and dimM
for m in range(0, 6):
    A = rng.standard_normal((2, 2))
    Ap = [[const(A[r, s]) for s in range(2)] for r in range(2)]
    tot = sum(dimM(m, j2) * np.trace(np.array([[e[0] for e in row] for row in pi_j(Ap, m, j2)]))
              for j2 in range(m % 2, m + 1, 2))
    ok &= abs(tot - np.trace(A) ** m) < 1e-8 * max(1, abs(np.trace(A)) ** m)
print('character identity:', ok)

for n in (1, 2, 3):
    D = 3 ** n
    Qm, _ = np.linalg.qr(rand_c((D, 2), rng))
    Pi = Qm @ Qm.conj().T
    mom = moments(Pi, n)
    c = coefficients(mom, n)
    worst_psd, worst_rank2, worst_comm = 0.0, 0.0, 0.0
    tot = 0
    for lamv in itertools.product('tsa', repeat=n):
        B = block(c, lamv)
        ev = np.linalg.eigvalsh((B + B.conj().T) / 2)
        worst_psd = min(worst_psd, ev.min())
        bcount = sum(1 for l in lamv if l == 'a')
        # sign projector on this block: (1/6) sum sgn(p) (x) rho_lam(p)
        P = np.zeros_like(B)
        for p in PERMS:
            M = np.array([[1.0]])
            for l in lamv:
                M = np.kron(M, rho(l, p))
            P += sign(p) * M / 6
        worst_rank2 = max(worst_rank2, np.abs(B @ P).max())
        for p in PERMS:
            M = np.array([[1.0]])
            for l in lamv:
                M = np.kron(M, rho(l, p))
            worst_comm = max(worst_comm, np.abs(B @ M - M @ B).max())
        tot += np.prod([DIMV[l] for l in lamv]) * np.trace(B).real
    print(f'n={n}: min eig of blocks {worst_psd:.2e}, |T P_sign| {worst_rank2:.1e}, |[T, tau]| {worst_comm:.1e}, '
          f'Tr T {tot:.6f} (expect 8)')
    ok &= worst_psd > -1e-9 and worst_rank2 < 1e-9 and worst_comm < 1e-9 and abs(tot - 8) < 1e-8
    # functionals: objective and A_j, B_j through blocks vs direct
    def via_blocks(f):          # f: dict perm -> float
        val = 0
        for lamv in itertools.product('tsa', repeat=n):
            B = block(c, lamv)
            F = np.array([[1.0]])
            for l in lamv:
                F = np.kron(F, sum(cf * rho(l, p) for p, cf in f.items()))
            val += np.prod([DIMV[l] for l in lamv]) * np.trace(B @ F)
        return val.real
    obj = via_blocks({(1, 0, 2): 2.0, (0, 1, 2): -1.0}) / 2 / 2 ** n
    ok &= abs(obj - q_dense(Pi, n).real) < 1e-9
    print(f'   q_n via blocks {obj:.10f}, direct {q_dense(Pi, n).real:.10f}')
print('PASS' if ok else 'FAIL')
