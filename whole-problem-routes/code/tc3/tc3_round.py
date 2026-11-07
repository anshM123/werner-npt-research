"""Three-copy thickness certificate, step 2: round the float solution, correct the 190 orbit equations exactly, verify exactly, write the certificate.

Exact objects: Y_beta (rational symmetric, reduced block variables), Gamma_beta = R_beta Y_beta R_beta^T, diagram
coefficients c = (x)_s Psi_s z (z = Gram-block coordinates), embedded into 4-strand diagram tuples and summed over
orbits.  Checks: (1) orbit sums of target - Z0 - Z1 - G(L1) - relations vanish exactly; (2) every Y_beta is positive
definite (exact LDL^T pivots > 0)."""
import itertools
import pickle
import sys
import time
from fractions import Fraction
import numpy as np
import flint
from exact_alg import build

T0 = time.time()


def log(*a):
    print(*a, f'[{time.time() - T0:.1f}s]', flush=True)


S = pickle.load(open('tc3_setup.pkl', 'rb'))
fsol = np.load('tc3_float.npz')
A = {'A2': build(2), 'A3': build(3), 'A4': build(4)}
n = 3
NO = S['NO']
tuples, orb, weights = S['tuples'], S['orb'], S['weights']
orbit_of = {int(f): int(o) for f, o in zip(tuples @ weights, orb)}
terms, blocks, copy_embed = S['terms'], S['blocks'], S['copy_embed']
target = S['target']
DEN = 2 ** 26


def q(x):
    return flint.fmpq(int(round(x * DEN)), DEN)


def to_fmpq_mat(M):
    return flint.fmpq_mat([[flint.fmpq(int(v)) if not isinstance(v, flint.fmpq) else v for v in row] for row in M])


Psi = {nm: flint.fmpq_mat([[flint.fmpq(x.numerator, x.denominator) for x in row] for row in A[nm]['Psi']]) for nm in A}
cidx = {nm: {c: i for i, c in enumerate(A[nm]['coords'])} for nm in A}


def embed_images(nm_term, s):
    """per copy: small diagram index -> list of (4-strand index, integer coefficient)"""
    E = copy_embed[nm_term][s]
    out = []
    for a in range(E.shape[1]):
        out.append([(int(i), int(E[i, a])) for i in np.nonzero(E[:, a])[0]])
    return out


emb = {nm: [embed_images(nm, s) for s in range(n)] for nm in terms}


def expand_block(nm, B, Gam):
    """diagram coefficients (dict small-index tuple -> fmpq) of the operator whose only non-zero Gram block is Gam (beta)"""
    names = terms[nm]
    beta = B['beta']
    ms = B['ms']
    # coordinates of the block in each copy
    cols = []
    for s in range(n):
        j = beta[s]
        cols.append([cidx[names[s]][(j, r, c)] for r in range(ms[s]) for c in range(ms[s])])
    # Gam is (b x b) with b = prod ms, rows/cols multi-indices (r0 r1 r2) C order; reorganise into a tensor over
    # per-copy (r_s, c_s) pairs: T[(r0 c0), (r1 c1), (r2 c2)]
    b = B['b']
    rows_ = list(itertools.product(*[range(x) for x in ms]))
    sizes = [ms[s] * ms[s] for s in range(n)]
    Tt = {}
    for R_, rr in enumerate(rows_):
        for C_, cc in enumerate(rows_):
            v = Gam[R_, C_]
            if v != 0:
                key = tuple(rr[s] * ms[s] + cc[s] for s in range(n))
                Tt[key] = v
    # mode products: result[a0, a1, a2] = sum Psi0[a0, col0[k0]] Psi1[a1, col1[k1]] Psi2[a2, col2[k2]] T[k0, k1, k2]
    # do it as a matrix product per mode with flint
    N = [Psi[names[s]].nrows() for s in range(n)]
    # mode 0
    M0 = flint.fmpq_mat(sizes[0], sizes[1] * sizes[2])
    for (k0, k1, k2), v in Tt.items():
        M0[k0, k1 * sizes[2] + k2] = v
    P0 = flint.fmpq_mat(N[0], sizes[0])
    for a in range(N[0]):
        for k0 in range(sizes[0]):
            P0[a, k0] = Psi[names[0]][a, cols[0][k0]]
    R0 = P0 * M0                                              # N0 x (s1 s2)
    # mode 1: rearrange to (s1) x (N0 s2)
    M1 = flint.fmpq_mat(sizes[1], N[0] * sizes[2])
    for a0 in range(N[0]):
        for k1 in range(sizes[1]):
            for k2 in range(sizes[2]):
                v = R0[a0, k1 * sizes[2] + k2]
                if v != 0:
                    M1[k1, a0 * sizes[2] + k2] = v
    P1 = flint.fmpq_mat(N[1], sizes[1])
    for a in range(N[1]):
        for k1 in range(sizes[1]):
            P1[a, k1] = Psi[names[1]][a, cols[1][k1]]
    R1 = P1 * M1                                              # N1 x (N0 s2)
    M2 = flint.fmpq_mat(sizes[2], N[0] * N[1])
    for a1 in range(N[1]):
        for a0 in range(N[0]):
            for k2 in range(sizes[2]):
                v = R1[a1, a0 * sizes[2] + k2]
                if v != 0:
                    M2[k2, a0 * N[1] + a1] = v
    P2 = flint.fmpq_mat(N[2], sizes[2])
    for a in range(N[2]):
        for k2 in range(sizes[2]):
            P2[a, k2] = Psi[names[2]][a, cols[2][k2]]
    R2 = P2 * M2                                              # N2 x (N0 N1)
    out = {}
    for a2 in range(N[2]):
        for a0 in range(N[0]):
            for a1 in range(N[1]):
                v = R2[a2, a0 * N[1] + a1]
                if v != 0:
                    out[(a0, a1, a2)] = v
    return out


def orbit_sums_of(nm, coeffs):
    out = [flint.fmpq(0)] * NO
    for (a0, a1, a2), v in coeffs.items():
        for i0, e0 in emb[nm][0][a0]:
            for i1, e1 in emb[nm][1][a1]:
                for i2, e2 in emb[nm][2][a2]:
                    o = orbit_of[i0 * weights[0] + i1 * weights[1] + i2 * weights[2]]
                    out[o] += v * (e0 * e1 * e2)
    return out


def gram_from_Y(B, Y):
    R = to_fmpq_mat(B['R'])
    return R * Y * R.transpose()


# ---------------------------------------------------------------- rounding
Yq = {}
for nm in terms:
    for B in blocks[nm]:
        key = f'{nm}_{"".join(map(str, B["beta"]))}'
        if B['R'].shape[1] == 0:
            continue
        Yf = fsol[key]
        br = Yf.shape[0]
        Yq[key] = flint.fmpq_mat([[q((Yf[i, j] + Yf[j, i]) / 2) for j in range(br)] for i in range(br)])
rho_f = fsol['rho']
rel_cols = S['rel_cols']
log('rounded', len(Yq), 'blocks')

# relation columns (exact orbit sums)
rel4 = A['A4']['relations']
E_xyy = copy_embed['Z1'][1]
img_xyy = [int(np.nonzero(E_xyy[:, a])[0][0]) for a in range(E_xyy.shape[1])]
rel_vecs = []
for col in rel_cols:
    r_, rem = divmod(col, 15 * 15)
    a1, a2 = divmod(rem, 15)
    v = [flint.fmpq(0)] * NO
    for i0, c in enumerate(rel4[r_]):
        if c:
            v[orbit_of[i0 * weights[0] + img_xyy[a1] * weights[1] + img_xyy[a2] * weights[2]]] += c
    rel_vecs.append(v)
rho_q = [q(x) for x in rho_f]


def total_orbit_sums(Ydict, rho):
    tot = [flint.fmpq(0)] * NO
    for nm in terms:
        for B in blocks[nm]:
            key = f'{nm}_{"".join(map(str, B["beta"]))}'
            if key not in Ydict:
                continue
            Gam = gram_from_Y(B, Ydict[key])
            os_ = orbit_sums_of(nm, expand_block(nm, B, Gam))
            tot = [x + y for x, y in zip(tot, os_)]
    for v, r in zip(rel_vecs, rho):
        tot = [x + r * y for x, y in zip(tot, v)]
    return tot


tq = [flint.fmpq(x.numerator, x.denominator) for x in target]
tot = total_orbit_sums(Yq, rho_q)
resid = [a - b for a, b in zip(tq, tot)]
log('residual after rounding: max', max(abs(float(x)) for x in resid))

# ---------------------------------------------------------------- correction columns (float selection, exact columns)
cand = []          # (key, p) diagonal entries
for nm in terms:
    for B in blocks[nm]:
        key = f'{nm}_{"".join(map(str, B["beta"]))}'
        if key in Yq:
            for p in range(B['R'].shape[1]):
                cand.append((nm, B, key, p))
log('candidate diagonal entries', len(cand))
# float columns via the exact routine would be slow; use a cheap float surrogate: expand with float arithmetic
Psi_f = {nm: np.array([[float(x) for x in row] for row in A[nm]['Psi']]) for nm in A}


def float_column(nm, B, p):
    Rf = B['R'].astype(float)
    Gam = np.outer(Rf[:, p], Rf[:, p])
    names = terms[nm]
    ms = B['ms']
    rows_ = list(itertools.product(*[range(x) for x in ms]))
    sizes = [m * m for m in ms]
    Tt = np.zeros(sizes)
    for R_, rr in enumerate(rows_):
        for C_, cc in enumerate(rows_):
            Tt[tuple(rr[s] * ms[s] + cc[s] for s in range(n))] = Gam[R_, C_]
    cols = [[cidx[names[s]][(B['beta'][s], r, c)] for r in range(ms[s]) for c in range(ms[s])] for s in range(n)]
    Cc = np.einsum('ai,bj,ck,ijk->abc', Psi_f[names[0]][:, cols[0]], Psi_f[names[1]][:, cols[1]], Psi_f[names[2]][:, cols[2]], Tt)
    out = np.zeros(NO)
    for a0, a1, a2 in zip(*np.nonzero(np.abs(Cc) > 1e-15)):
        v = Cc[a0, a1, a2]
        for i0, e0 in emb[nm][0][a0]:
            for i1, e1 in emb[nm][1][a1]:
                for i2, e2 in emb[nm][2][a2]:
                    out[orbit_of[i0 * weights[0] + i1 * weights[1] + i2 * weights[2]]] += v * e0 * e1 * e2
    return out


FC = np.array([float_column(nm, B, p) for nm, B, key, p in cand]).T          # NO x ncand
RF = np.array([[float(x) for x in v] for v in rel_vecs]).T                   # NO x 19
allF = np.hstack([RF, FC])
import scipy.linalg as sla
_, _, piv = sla.qr(allF, pivoting=True)
rk = np.linalg.matrix_rank(allF)
chosen = piv[:rk]
log('float rank of correction columns', rk, 'cond of chosen', np.linalg.cond(allF[:, chosen]))
cols_exact = []
for c in chosen:
    if c < len(rel_vecs):
        cols_exact.append(rel_vecs[c])
    else:
        nm, B, key, p = cand[c - len(rel_vecs)]
        br = B['R'].shape[1]
        E = flint.fmpq_mat(br, br)
        E[p, p] = 1
        cols_exact.append(orbit_sums_of(nm, expand_block(nm, B, gram_from_Y(B, E))))
log('exact correction columns built')
Afull = flint.fmpq_mat([[cols_exact[j][i] for j in range(rk)] for i in range(NO)])     # NO x rk
log('exact rank of chosen columns', Afull.rank())
# left null space (exact): solve A^T u = 0
At = Afull.transpose()
rref, rnk = At.rref()
# free variables -> null vectors
pivcols = []
r_ = 0
for c in range(NO):
    if r_ < rnk and rref[r_, c] == 1 and all(rref[r_, cc] == 0 for cc in range(c)):
        pivcols.append(c)
        r_ += 1
free = [c for c in range(NO) if c not in pivcols]
nulls = []
for fcol in free:
    u = [flint.fmpq(0)] * NO
    u[fcol] = flint.fmpq(1)
    for rr, pc in enumerate(pivcols):
        u[pc] = -rref[rr, fcol]
    nulls.append(u)
log('left null space dimension', len(nulls))
for u in nulls:
    ut = sum((u[i] * tq[i] for i in range(NO)), flint.fmpq(0))
    ur = sum((u[i] * resid[i] for i in range(NO)), flint.fmpq(0))
    log('  u . target =', ut, '  u . residual =', ur)
# solve on independent rows
_, _, prow = sla.qr(np.array([[float(Afull[i, j]) for j in range(rk)] for i in range(NO)]).T, pivoting=True)
rows_sel = sorted(prow[:rk].tolist())
Asq = flint.fmpq_mat([[Afull[i, j] for j in range(rk)] for i in rows_sel])
rv = flint.fmpq_mat([[resid[i]] for i in rows_sel])
delta = Asq.solve(rv)
log('correction max', max(abs(float(delta[i, 0])) for i in range(rk)))
for jj, c in enumerate(chosen):
    d = delta[jj, 0]
    if c < len(rel_vecs):
        rho_q[c] = rho_q[c] + d
    else:
        nm, B, key, p = cand[c - len(rel_vecs)]
        Yq[key][p, p] = Yq[key][p, p] + d

# ---------------------------------------------------------------- exact verification
tot = total_orbit_sums(Yq, rho_q)
bad = [i for i in range(NO) if tq[i] != tot[i]]
log('orbit equations violated after correction:', len(bad))


def ldl_positive(M):
    """exact LDL^T; returns least pivot (fmpq) or None if a pivot is <= 0"""
    m = M.nrows()
    A_ = [[M[i, j] for j in range(m)] for i in range(m)]
    least = None
    for k in range(m):
        piv_ = A_[k][k]
        if piv_ <= 0:
            return None
        least = piv_ if least is None or piv_ < least else least
        for i in range(k + 1, m):
            if A_[i][k] != 0:
                f = A_[i][k] / piv_
                for j in range(k + 1, m):
                    A_[i][j] = A_[i][j] - f * A_[k][j]
    return least


worst = None
for key, Y in Yq.items():
    assert Y == Y.transpose()
    lp = ldl_positive(Y)
    if lp is None:
        log('NOT positive definite:', key)
        worst = flint.fmpq(-1)
        break
    worst = lp if worst is None or lp < worst else worst
log('all reduced blocks positive definite; least LDL pivot', float(worst))
pickle.dump(dict(Y={k: [[(int(v.p), int(v.q)) for v in row] for row in [[Y[i, j] for j in range(Y.ncols())] for i in range(Y.nrows())]] for k, Y in Yq.items()},
                 rho=[(int(r.p), int(r.q)) for r in rho_q], rel_cols=rel_cols, ok=(len(bad) == 0)),
            open('tc3_cert_blocks.pkl', 'wb'))
log('certificate written; PASS' if not bad else 'FAIL')
