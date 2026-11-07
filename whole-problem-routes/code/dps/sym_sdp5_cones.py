"""Level-2 DPS certificate SDP with the full copy-permutation (S_n) and transpose symmetry, for general n.

Same certificate format as level2_sdp.py:
    P (W (x) I_B' - t I) P = Z + Gamma_A(Q_A) + P Gamma_B(Q_B) P,   Q_A, Q_B, Z >= 0, Q_A = P Q_A P, Z = P Z P.
Per-copy block bases U_j are chosen inside the S = +1 / S = -1 eigenspaces of SWAP(B,B') (S diagonal), so the
symmetric range of every block is a set of coordinates.  Unknowns: the block-coordinate tensors of Q_A, Q_B,
which are invariant under permutations of the copies and under transposition; one variable per orbit.
Only blocks beta with j_1 <= ... <= j_n are imposed PSD (the others are permutation conjugates).

usage: python sym_sdp.py n a_num a_den [mode] [tfix] [eps]
   mode 'maxt' (default) or 'center' (t = tfix fixed, maximise the least eigenvalue s of all blocks)
"""
import sys
import time
import itertools
import numpy as np
import cvxpy as cp
import brauer as br

n = int(sys.argv[1])
anum, aden = int(sys.argv[2]), int(sys.argv[3])
mode = sys.argv[4] if len(sys.argv) > 4 else 'maxt'
tfix = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
eps = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-9
a = anum / aden
T0 = time.time()
SP = [0, 1, 2, 3]
NB = br.NDIAG


def log(*x):
    print(*x, f'[{time.time() - T0:.1f}s]', flush=True)


# ---------------------------------------------------------------- per copy
Bex = br.spin_m0_bases()
Dm = [M.astype(float) for M in br.DMAT]
nm = br.named()
S = Dm[nm['S_BB']]
U, sgn = {}, {}
for j in SP:
    Bj = np.array(Bex[j].tolist(), dtype=float)
    parts, sg = [], []
    for sign in (+1, -1):
        X = (Bj + sign * S @ Bj) / 2
        u_, s_, vt = np.linalg.svd(X, full_matrices=False)
        r = int((s_ > 1e-9).sum())
        if r:
            parts.append(u_[:, :r])
            sg += [sign] * r
    U[j] = np.hstack(parts)
    sgn[j] = np.array(sg)
    assert U[j].shape[1] == Bj.shape[1]
    assert np.allclose(U[j].T @ S @ U[j], np.diag(sgn[j]))
m = {j: U[j].shape[1] for j in SP}
coords = []          # per-copy block coordinates: (j, r, s)
for j in SP:
    for r in range(m[j]):
        for s in range(m[j]):
            coords.append((j, r, s))
assert len(coords) == NB
cindex = {c: i for i, c in enumerate(coords)}
tau = np.array([cindex[(j, s, r)] for (j, r, s) in coords])        # transpose of the block entry
Phi = np.zeros((NB, NB))
for a_ in range(NB):
    for i, (j, r, s) in enumerate(coords):
        Phi[i, a_] = U[j][:, r] @ Dm[a_] @ U[j][:, s]
Phinv = np.linalg.inv(Phi)


def ptmat(strands):
    pi = br.pt_perm(strands)
    T = np.zeros((NB, NB))
    for a_ in range(NB):
        T[pi[a_], a_] = 1.0
    return Phi @ T @ Phinv


MA, MB = ptmat([0]), ptmat([1])
cdiag = np.zeros(NB)
cdiag[nm['I']] = a
cdiag[nm['F_AB']] -= 0.25
cdiag[nm['Om_AB']] -= 0.25
idd = np.zeros(NB)
idd[nm['I']] = 1.0
cblk, iblk = Phi @ cdiag, Phi @ idd
jof = np.array([c[0] for c in coords])
rsg = np.array([sgn[c[0]][c[1]] for c in coords])      # S-sign of the row index
csg = np.array([sgn[c[0]][c[2]] for c in coords])      # S-sign of the column index

# ---------------------------------------------------------------- global coordinates and orbits
allk = np.array(list(itertools.product(range(NB), repeat=n)), dtype=np.int64)    # (15^n, n), C order
N = allk.shape[0]
pw = NB ** np.arange(n - 1, -1, -1)
canon = np.minimum(np.sort(allk, axis=1) @ pw, np.sort(tau[allk], axis=1) @ pw)   # S_n x transpose orbits
orbit_reps, orbit_of = np.unique(canon, return_inverse=True)
NO = len(orbit_reps)
symsym = (np.prod(rsg[allk], axis=1) == 1) & (np.prod(csg[allk], axis=1) == 1)
log(f'n={n} a={a}: {N} coordinates, {NO} orbits; Sym x Sym coordinates {symsym.sum()}')
orbA = np.unique(orbit_of[symsym])                       # orbits used by Q_A and Z
posA = -np.ones(NO, dtype=np.int64)
posA[orbA] = np.arange(len(orbA))
NA = len(orbA)


def kron_row(M, k):
    """row of M^{(x)n} at global coordinate k (length 15^n)"""
    r = np.array([1.0])
    for s in range(n):
        r = np.kron(r, M[k[s]])
    return r


repk = allk[orbit_reps]                                   # representative tuple of each orbit (the canonical one)
# Z at Sym x Sym orbits: z = c^{(x)n} - t i^{(x)n} - MA^{(x)n} yA - MB^{(x)n} yB, yA = E_A uA, yB = E_B uB
RA = np.zeros((NA, NA))
RB = np.zeros((NA, NO))
cz = np.zeros(NA)
iz = np.zeros(NA)
for ii, o in enumerate(orbA):
    k = repk[o]
    rowA = kron_row(MA, k)
    rowB = kron_row(MB, k)
    RB[ii] = np.bincount(orbit_of, weights=rowB, minlength=NO)
    fullA = np.bincount(orbit_of, weights=rowA * symsym, minlength=NO)
    RA[ii] = fullA[orbA]
    cz[ii] = np.prod(cblk[k])
    iz[ii] = np.prod(iblk[k])
log('dense maps built', RA.shape, RB.shape)

# ---------------------------------------------------------------- blocks
blocks = []
for beta in itertools.combinations_with_replacement(SP, n):
    ms = [m[j] for j in beta]
    rows = list(itertools.product(*[range(x) for x in ms]))
    b = len(rows)
    flat = np.zeros((b, b), dtype=np.int64)
    for R, rr in enumerate(rows):
        for C_, cc in enumerate(rows):
            f = 0
            for s in range(n):
                f = f * NB + cindex[(beta[s], rr[s], cc[s])]
            flat[R, C_] = f
    rs = np.array([np.prod([sgn[beta[s]][rr[s]] for s in range(n)]) for rr in rows])
    sym = np.where(rs == 1)[0]
    blocks.append(dict(beta=beta, b=b, flat=flat, sym=sym))
log('blocks', len(blocks), 'max b', max(B['b'] for B in blocks), 'sum b^2', sum(B['b'] ** 2 for B in blocks))

uA = cp.Variable(NA)
uB = cp.Variable(NO)
t = cp.Variable() if mode == 'maxt' else tfix
s_ = cp.Variable() if mode == 'center' else 0.0
zv = cp.Variable(NA)
cons = [zv == cz - t * iz - RA @ uA - RB @ uB]
z = zv
for B in blocks:
    ob = orbit_of[B['flat']]
    Xb = cp.reshape(uB[ob.reshape(-1)], (B['b'], B['b']), order='C')
    cons.append((Xb + Xb.T) / 2 >> s_ * np.eye(B['b']))
    if len(B['sym']):
        sub = ob[np.ix_(B['sym'], B['sym'])]
        k = len(B['sym'])
        Xa = cp.reshape(uA[posA[sub].reshape(-1)], (k, k), order='C')
        Xz = cp.reshape(z[posA[sub].reshape(-1)], (k, k), order='C')
        assert (posA[sub] >= 0).all()
        cons.append((Xa + Xa.T) / 2 >> s_ * np.eye(k))
        cons.append((Xz + Xz.T) / 2 >> s_ * np.eye(k))
cones = sys.argv[8] if len(sys.argv) > 8 else 'AB'
if 'A' not in cones:
    cons.append(uA == 0)
if 'B' not in cones:
    cons.append(uB == 0)
obj = cp.Maximize(t if mode == 'maxt' else s_)
prob = cp.Problem(obj, cons)
log('problem built; variables', 2 * NA + NO)
solver = sys.argv[7] if len(sys.argv) > 7 else 'CLARABEL'
if solver == 'CLARABEL':
    prob.solve(solver=cp.CLARABEL, verbose=True, max_iter=400, tol_gap_abs=eps, tol_gap_rel=eps, tol_feas=eps)
else:
    prob.solve(solver=cp.SCS, verbose=True, max_iters=100000, eps_abs=eps, eps_rel=eps)
log('status', prob.status, 'value', prob.value)
uAv, uBv = uA.value, uB.value
# expand to full block-coordinate tensors and then to diagram coordinates
yA = np.zeros(NO)
yA[orbA] = uAv
yA = yA[orbit_of] * symsym
yB = uBv[orbit_of]


def to_diagram(y):
    T = y.reshape([NB] * n)
    for s in range(n):
        T = np.moveaxis(np.tensordot(Phinv, T, axes=([1], [s])), 0, s)
    return T.reshape(-1)


xA, xB = to_diagram(yA), to_diagram(yB)
tag = f'sym_n{n}_{anum}_{aden}_{mode}_' + (sys.argv[8] if len(sys.argv) > 8 else 'AB')
np.savez(f'sol_{tag}.npz', xA=xA, xB=xB, t=(t.value if mode == 'maxt' else tfix),
         s=(s_.value if mode == 'center' else 0.0))
log('saved', f'sol_{tag}.npz', 't =', (t.value if mode == 'maxt' else tfix), 's =',
    (s_.value if mode == 'center' else None))
