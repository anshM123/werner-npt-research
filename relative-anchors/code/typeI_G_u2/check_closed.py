"""Exact checks of the closed forms of closed_u2.pkl (u = 2) against direct exact elimination of the whole type I
block (cf.TypeI over Fractions), for m = 1, 2, 3, 4, 5 and random rational weights (so the z-scaling is rational).
G_tt is compared with G00(s_t, sigma - s_t, 0, s_3, s_4) and G_tt' with G01(s_t, s_t', sigma - s_t - s_t', s_3, s_4):
this is the pairwise structure (S4) at these points.  Also checks, in floating point, that tr(Omega_t G Omega_t' G^T)
equals ||Ghat_tt'||_F^2 and that lambda_min(Ghat_tt) is the root used in (B), with Ghat built as in sweep.py of the
repo (orthonormal frame, aligned direction rescaled by c = sqrt(1 + S/sigma))."""
import sys, os, random, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import numpy as np
import flint
from cf import TypeI

ctx5 = flint.fmpq_mpoly_ctx.get(('s0', 's1', 's2', 's3', 's4'), 'lex')
D = pickle.load(open('closed_u2.pkl', 'rb'))
P = {nm: [[(ctx5.from_dict(D[nm][i][j][0]), ctx5.from_dict(D[nm][i][j][1])) for j in range(3)] for i in range(3)] for nm in D}


def closed(nm, s):
    q = [flint.fmpq(x.numerator, x.denominator) for x in s]
    out = []
    for i in range(3):
        row = []
        for j in range(3):
            n_, d_ = P[nm][i][j]
            v = n_(*q) / d_(*q)
            row.append(Fr(int(v.p), int(v.q)))
        out.append(row)
    return out


def schur_exact(w, T):
    ti = TypeI(w, T, F=Fr)
    M = ti.build()
    U = ti.U; r0 = U[0]
    tb = []; scale = []
    for t in ti.T:
        tb.append(ti.idx[('R', t)]); scale.append(Fr(1))
        for h in U:
            tb.append(ti.idx[('Z', t, h)]); scale.append(w[t] * w[h])
    el = [i for i in range(ti.n) if i not in tb and not (ti.labels[i][0] == 'P' and ti.labels[i][1] == r0)]
    order = el + tb
    A = [[M[i][j] for j in order] for i in order]
    ne = len(el)
    for p in range(ne):
        pv = A[p][p]
        assert pv > 0
        for i in range(p + 1, len(order)):
            if A[i][p] == 0:
                continue
            f = A[i][p] / pv
            for j in range(p + 1, len(order)):
                A[i][j] -= f * A[p][j]
    nt = len(tb)
    return [[A[ne + i][ne + j] * scale[i] * scale[j] for j in range(nt)] for i in range(nt)], ti


random.seed(7)
ok = True
ncmp = 0
for m in (1, 2, 3, 4, 5):
    for trial in range(3 if m < 5 else 1):
        wT = [Fr(random.randint(1, 30), random.randint(1, 30)) for _ in range(m)]
        wU = [Fr(random.randint(1, 30), random.randint(1, 30)) for _ in range(2)]
        w = wT + wU
        G, ti = schur_exact(w, list(range(m)))
        s = [x * x for x in w]
        sig = sum(s[:m]); s3, s4 = s[m], s[m + 1]
        for t in range(m):
            ref = closed('G00', [s[t], sig - s[t], Fr(0), s3, s4])
            blk = [[G[3 * t + i][3 * t + j] for j in range(3)] for i in range(3)]
            ok &= blk == ref; ncmp += 1
            for tp in range(m):
                if tp == t:
                    continue
                ref = closed('G01', [s[t], s[tp], sig - s[t] - s[tp], s3, s4])
                blk = [[G[3 * t + i][3 * tp + j] for j in range(3)] for i in range(3)]
                ok &= blk == ref; ncmp += 1
        print(f'm = {m}, weights {[str(x) for x in w]}: blocks equal to the closed forms: {ok}', flush=True)
print(f'{ncmp} blocks compared exactly')

# float cross-check of the norm formulas with the frame of sweep.py
rng = np.random.default_rng(3)
err = 0.0
for trial in range(20):
    m = 3
    w = np.exp(rng.uniform(-1.5, 1.5, m + 2))
    G, ti = schur_exact([Fr(x).limit_denominator(10 ** 6) for x in w], list(range(m)))
    wf = np.array([float(Fr(x).limit_denominator(10 ** 6)) for x in w]); s = wf ** 2
    Gz = np.array([[float(x) for x in row] for row in G])
    sig = s[:m].sum(); S = s[m:].sum(); wU = wf[m:]; wh = wU / np.linalg.norm(wU); perp = np.array([-wh[1], wh[0]])
    c = np.sqrt(1 + S / sig)
    def frame(t):
        Z = np.diag([1, wf[t] * wU[0], wf[t] * wU[1]])     # zeta = Z z
        C = np.zeros((3, 3)); C[0, 0] = 1; C[1:, 1] = c * wh; C[1:, 2] = perp   # zeta = C (rho, a, b)
        return np.linalg.solve(Z, C)                         # z = Z^-1 C (rho, a, b)
    for t in range(m):
        Ft = frame(t)
        Ght = Ft.T @ Gz[3 * t:3 * t + 3, 3 * t:3 * t + 3] @ Ft
        Om = np.diag([1.0, 0, 0]); Om[1:, 1:] = (np.diag(1 / s[m:]) + 1 / sig) / s[t]
        for tp in range(m):
            if tp == t:
                continue
            Gb = Gz[3 * t:3 * t + 3, 3 * tp:3 * tp + 3]
            Ftp = frame(tp)
            Ghat = Ft.T @ Gb @ Ftp
            Omp = np.diag([1.0, 0, 0]); Omp[1:, 1:] = (np.diag(1 / s[m:]) + 1 / sig) / s[tp]
            err = max(err, abs(np.linalg.norm(Ghat, 'fro') ** 2 - np.trace(Om @ Gb @ Omp @ Gb.T)) / max(1e-300, np.linalg.norm(Ghat, 'fro') ** 2))
        # (B): Ghat_tt - h I > 0  <=>  G_z - h Z D^-2 Z > 0 ; compare smallest generalized eigenvalue
        lam = np.linalg.eigvalsh(Ght)[0]
        Oinv = np.zeros((3, 3)); Oinv[0, 0] = 1; Oinv[1:, 1:] = s[t] * (np.diag(s[m:]) - np.outer(s[m:], s[m:]) / (sig + S))
        Gt = Gz[3 * t:3 * t + 3, 3 * t:3 * t + 3]
        L = np.linalg.cholesky(Oinv); Li = np.linalg.inv(L)
        lam2 = np.linalg.eigvalsh(Li @ Gt @ Li.T)[0]
        err = max(err, abs(lam - lam2) / abs(lam))
print(f'float: max relative mismatch of the norm/eigenvalue formulas {err:.2e}')
ok &= err < 1e-8
print('PASS' if ok else 'FAIL')
