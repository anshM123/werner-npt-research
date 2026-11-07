"""Level-2 DPS for block-positivity of C_{a_1} (x) ... (x) C_{a_n} with one parameter per copy (as level2_sdp.py, which
has a single parameter).  Floating point (cvxpy).  Used for the vertex (1/2, 3/4, 3/4) of Conjecture QT_m (Proposition 7.20),
where the best t is negative (about -4.7e-3), so this level does not certify 2-positivity of P_{3/4} (x) P_{3/4}.

usage: python level2_multi.py a_1,a_2,...,a_n [mode] [tfix] [solver] [eps] [max_iters]
   mode 'maxt' (default) or 'center';  solver CLARABEL (default linear solver) or SCS (eps, max_iters for SCS)
   e.g.  python level2_multi.py 0.5,0.75,0.75 maxt 0 SCS 1e-5 60000      (about 25 minutes)
         python level2_multi.py 0.5,0.75 maxt 0 CLARABEL                  (control: t* = 0 at the two-copy vertex)
"""
import sys
import time
import numpy as np
import scipy.sparse as sps
import cvxpy as cp
import brauer as br
import blocks as bl

alist = [float(v) for v in sys.argv[1].split(',')]
n = len(alist)
mode = sys.argv[2] if len(sys.argv) > 2 else 'maxt'
tfix = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
a = alist[0]; anum = 0; aden = 1
t0 = time.time()
pc = bl.per_copy_float()
nm = pc['named']
N = br.NDIAG ** n

# per-copy C_a (x) I_B' in diagram coordinates
idd = np.zeros(br.NDIAG)
idd[nm['I']] = 1.0
iblk = pc['Phi'] @ idd
c3 = np.array([1.0])
for aa in alist:
    cdiag = np.zeros(br.NDIAG)
    cdiag[nm['I']] = aa
    cdiag[nm['F_AB']] -= 0.25
    cdiag[nm['Om_AB']] -= 0.25
    c3 = np.kron(c3, pc['Phi'] @ cdiag)
i3 = bl.coeff_power(iblk, n)

MA = pc['ptmat']([0])
MB = pc['ptmat']([1])
KA = bl.kron_chain_mats(MA, n)
KB = bl.kron_chain_mats(MB, n)
blks = bl.global_blocks(pc, n)
print('n', n, 'a', alist, 'blocks', len(blks), 'sum b^2', sum(B['b'] ** 2 for B in blks),
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
solver = sys.argv[4] if len(sys.argv) > 4 else 'CLARABEL'
if solver == 'CLARABEL':
    prob.solve(solver=cp.CLARABEL, verbose=True, max_iter=500, tol_gap_abs=1e-9, tol_gap_rel=1e-9, tol_feas=1e-9)
else:
    prob.solve(solver=cp.SCS, verbose=True, max_iters=int(sys.argv[6]) if len(sys.argv) > 6 else 20000, eps_abs=float(sys.argv[5]) if len(sys.argv) > 5 else 1e-7, eps_rel=float(sys.argv[5]) if len(sys.argv) > 5 else 1e-7)
print('status', prob.status, 'value', prob.value, 'time', time.time() - t0, flush=True)
if mode == 'maxt':
    print('t* =', t.value, ' (aligned product value prod (a_s - 1/2) =', float(np.prod([x - 0.5 for x in alist])), ')')
# save solution in block coordinates
yAv = np.asarray(yA.value)
yBv = np.asarray(yB.value)
tag = 'v_' + '_'.join(str(x) for x in alist) + '_' + mode
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
