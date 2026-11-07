"""Level-2 DPS (symmetric extension on the output side, with partial transposes) for block-positivity of
W = C_a^{(x)n} across (A_1..A_n) | (B_1..B_n).  Floating point (cvxpy).

Find t max and PSD Q_A (supported on A (x) Sym(B,B')), Q_B, Z (supported on A (x) Sym(B,B')) with
    P (W (x) I_B' - t I) P = Z + Gamma_A(Q_A) + P Gamma_B(Q_B) P ,     P = I_A (x) P_sym(B,B').
All operators are taken O(3)^n invariant, so everything is done in block coordinates (blocks.py).

usage: python level2_sdp.py n a_num a_den [mode] [tfix]
   mode = 'maxt' (default) : maximise t
   mode = 'center'         : t fixed to tfix, maximise s with Q_A, Q_B, Z >= s I on their blocks
"""
import sys
import time
import numpy as np
import scipy.sparse as sps
import cvxpy as cp
import brauer as br
import blocks as bl

n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
anum = int(sys.argv[2]) if len(sys.argv) > 2 else 3
aden = int(sys.argv[3]) if len(sys.argv) > 3 else 4
mode = sys.argv[4] if len(sys.argv) > 4 else 'maxt'
tfix = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
a = anum / aden
t0 = time.time()
pc = bl.per_copy_float()
nm = pc['named']
N = br.NDIAG ** n

# per-copy C_a (x) I_B' in diagram coordinates
cdiag = np.zeros(br.NDIAG)
cdiag[nm['I']] = a
cdiag[nm['F_AB']] -= 0.25
cdiag[nm['Om_AB']] -= 0.25
idd = np.zeros(br.NDIAG)
idd[nm['I']] = 1.0
cblk = pc['Phi'] @ cdiag
iblk = pc['Phi'] @ idd
c3 = bl.coeff_power(cblk, n)
i3 = bl.coeff_power(iblk, n)

MA = pc['ptmat']([0])
MB = pc['ptmat']([1])
KA = bl.kron_chain_mats(MA, n)
KB = bl.kron_chain_mats(MB, n)
blks = bl.global_blocks(pc, n)
print('n', n, 'a', a, 'blocks', len(blks), 'sum b^2', sum(B['b'] ** 2 for B in blks),
      'sum bp^2', sum(B['bp'] ** 2 for B in blks), 'max b', max(B['b'] for B in blks), flush=True)

# selection matrices: global vector = sum_beta Sel_beta vec_C(Y_beta)
YB, YA, ZZ = [], [], []
cons = []
s = cp.Variable() if mode == 'center' else None
yB_terms, yA_terms = [], []
for B in blks:
    b, bp = B['b'], B['bp']
    rows = B['idx'].reshape(-1)
    Sel = sps.csr_matrix((np.ones(b * b), (rows, np.arange(b * b))), shape=(N, b * b))
    B['Sel'] = Sel
    Yb = cp.Variable((b, b), symmetric=True)
    YB.append(Yb)
    yB_terms.append(Sel @ cp.vec(Yb, order='C'))
    if bp > 0:
        Ya = cp.Variable((bp, bp), symmetric=True)
        Zz = cp.Variable((bp, bp), symmetric=True)
        YA.append(Ya)
        ZZ.append(Zz)
        Pi = B['Pi']
        yA_terms.append(Sel @ cp.vec(Pi @ Ya @ Pi.T, order='C'))
        if mode == 'center':
            cons += [Ya >> s * np.eye(bp), Zz >> s * np.eye(bp)]
        else:
            cons += [Ya >> 0, Zz >> 0]
    else:
        YA.append(None)
        ZZ.append(None)
    if mode == 'center':
        cons += [Yb >> s * np.eye(b)]
    else:
        cons += [Yb >> 0]
yB = cp.Variable(N)
yA = cp.Variable(N)
cons += [yB == sum(yB_terms), yA == sum(yA_terms)]


def chain(K, y):
    u = y
    for s_ in reversed(range(n)):
        v = cp.Variable(N)
        cons.append(v == K[s_] @ u)
        u = v
    return u


gA = chain(KA, yA)
gB = chain(KB, yB)
if mode == 'center':
    t = tfix
else:
    t = cp.Variable()
r = c3 - t * i3 - gA - gB
for B, Zz in zip(blks, ZZ):
    if Zz is None:
        continue
    b, Pi = B['b'], B['Pi']
    Rb = cp.reshape(r[B['idx'].reshape(-1)], (b, b), order='C')
    Ex = Zz - Pi.T @ ((Rb + Rb.T) / 2) @ Pi
    cons.append(cp.diag(Ex) == 0)
    if Zz.shape[0] > 1:
        cons.append(cp.upper_tri(Ex) == 0)
obj = cp.Maximize(s if mode == 'center' else t)
prob = cp.Problem(obj, cons)
print('built', time.time() - t0, flush=True)
solver = sys.argv[6] if len(sys.argv) > 6 else 'CLARABEL'
if solver == 'CLARABEL':
    prob.solve(solver=cp.CLARABEL, verbose=True, max_iter=500, tol_gap_abs=1e-9, tol_gap_rel=1e-9,
               tol_feas=1e-9, direct_solve_method='qdldl')
else:
    prob.solve(solver=cp.SCS, verbose=True, max_iters=20000, eps_abs=1e-7, eps_rel=1e-7)
print('status', prob.status, 'value', prob.value, 'time', time.time() - t0, flush=True)
if mode == 'maxt':
    print('t* =', t.value, ' (product-vector minimum (a-1/2)^n =', (a - 0.5) ** n, ')')
# save solution in block coordinates
yAv = np.asarray(yA.value)
yBv = np.asarray(yB.value)
tag = f'n{n}_{anum}_{aden}_{mode}'
np.savez(f'sol_{tag}.npz', yA=yAv, yB=yBv, t=(t.value if mode == 'maxt' else tfix),
         s=(s.value if mode == 'center' else 0.0))
# report min eigenvalues
mins = []
for B, Yb, Ya, Zz in zip(blks, YB, YA, ZZ):
    e = [np.linalg.eigvalsh(Yb.value).min()]
    if Ya is not None:
        e += [np.linalg.eigvalsh(Ya.value).min(), np.linalg.eigvalsh(Zz.value).min()]
    mins.append(min(e))
print('min eigenvalue over all blocks of Q_A, Q_B, Z:', min(mins))
