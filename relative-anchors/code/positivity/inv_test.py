"""Test 5: fixed-vector representations of matrices derived from the gauge-fixed block M_g:
  INV   W = M_g^{-1}                         (if W = sum c_e x_e x_e^T, Jacobi + Cauchy-Binet give SF principal minors of M_g)
  VAC   (vacuum block V of M_g)^{-1}, rescaled by diag(s^{N(y)}) (N(y) = letters 1 of the column y): Lemma G predicts
        a fixed-vector matrix I + (1/2) J^- C^{-1} J^+ in the coordinates eta = mu Psi
  SG1   Schur complement of M_g onto the states with a row letter (vacuum states eliminated)
  ROW a (principal submatrix of M_g on the states with row string a)^{-1}, rescaled by diag(s^{N(X)}) (all letters 1
        of the state X, row and column), the analogue of VAC for a row letter
For each: the exact obstruction of newton_obstr.py (with width certificate), the exact monomial-PSD test (q = 1), and
the LP search for a monomial rescaling followed by the exact monomial-PSD test.
Usage: python inv_test.py k [maxsize]"""
import sys, time
import flint
from psdtools import load_blocks
from flintfield import Field, F
from ylib import inv
from analyze import obstruction, mono_psd, rescale_lp, apply_rescale

k = int(sys.argv[1])
maxsize = int(sys.argv[2]) if len(sys.argv) > 2 else 99
t0 = time.time()


def toF(K, Mg):
    return [[F(K, x[0], x[1]) for x in row] for row in Mg]


def fromF(Mf):
    return [[(x.n, x.d) for x in row] for row in Mf]


def schur(Mf, keep, elim, K):
    A = [[Mf[a][b] for b in keep] for a in keep]
    if not elim:
        return A
    C = [[Mf[a][b] for b in elim] for a in elim]
    Ci = inv(C, K(0), K(1))
    Bm = [[Mf[a][b] for b in elim] for a in keep]
    out = []
    for i in range(len(keep)):
        row = []
        for j in range(len(keep)):
            v = A[i][j]
            for p in range(len(elim)):
                if Bm[i][p] == 0:
                    continue
                for q in range(len(elim)):
                    if Ci[p][q] != 0 and Bm[j][q] != 0:
                        v = v - Bm[i][p] * Ci[p][q] * Bm[j][q]
            row.append(v)
        out.append(row)
    return out


def report(tag, Mg, homog):
    ob = obstruction(Mg)
    if ob is None:
        obs = 'n/a (diagonal not subtraction-free)'
    else:
        obs = f'{len(ob[0])} obstructed pairs ({len(ob[1])} width)'
    nm, nf = mono_psd(Mg)
    d2 = rescale_lp(Mg, homog) if (ob is not None and ob[0] and not ob[1]) else None
    extra = ''
    if ob is not None and ob[1]:
        extra = '; no monomial rescaling can help (width certificate, exact)'
    elif ob is not None and ob[0]:
        if isinstance(d2, str):
            extra = '; LP skipped (too large)'
        elif d2 is None:
            extra = '; monomial-rescaling LP infeasible'
        else:
            nm2, nf2 = mono_psd(apply_rescale(Mg, d2, homog))
            extra = f'; LP rescaling 2*delta = {d2.tolist()}: non-PSD {nf2} of {nm2}'
    print(f'   {tag}: {obs}; monomial-PSD (q = 1): {nf} non-PSD of {nm}{extra}', flush=True)
    return ob, nf


homog = k in (2, 3)
for blk in load_blocks(k):
    n = len(blk['Mg'])
    if n > maxsize or n < 2:
        continue
    K = Field(list(blk['names']))
    Mf = toF(K, blk['Mg'])
    print(f'== {blk["name"]} ({blk["kind"]}), size {n}  [{time.time() - t0:.1f}s]', flush=True)
    W = inv(Mf, K(0), K(1))
    report('INV', fromF(W), homog)
    vac = [i for i, st in enumerate(blk['states']) if all(c == 0 for c in st[0])]
    rows = [i for i in range(n) if i not in vac]
    if len(vac) >= 2:
        V = [[Mf[a][b] for b in vac] for a in vac]
        Vi = inv(V, K(0), K(1))
        s = K.gens if not homog else K.gens
        if k == 4:
            s = [K(1)] + list(K.gens)
        dg = []
        for a in vac:
            y = blk['states'][a][1]
            m = K(1)
            for t, c in enumerate(y):
                if c == 1:
                    m = m * s[t]
            dg.append(m)
        Vr = [[Vi[i][j] * dg[i] * dg[j] for j in range(len(vac))] for i in range(len(vac))]
        report('VAC', fromF(Vr), homog)
    if vac and len(rows) >= 2:
        S = schur(Mf, rows, vac, K)
        report('SG1', fromF(S), homog)
    # row groups: states sharing one row string; inverse in the coordinates eta = mu Psi (rescaling by s^{N(X)})
    s_ = [K(1)] + list(K.gens) if k == 4 else list(K.gens)
    groups = {}
    for i, st in enumerate(blk['states']):
        if any(c != 0 for c in st[0]):
            groups.setdefault(st[0], []).append(i)
    for a_, G in sorted(groups.items()):
        if len(G) < 2:
            continue
        Vg = [[Mf[p][q] for q in G] for p in G]
        Vgi = inv(Vg, K(0), K(1))
        dg = []
        for p in G:
            a0, y0 = blk['states'][p]
            m = K(1)
            for t in range(k):
                for _ in range(int(a0[t] == 1) + int(y0[t] == 1)):
                    m = m * s_[t]
            dg.append(m)
        Vr = [[Vgi[i][j] * dg[i] * dg[j] for j in range(len(G))] for i in range(len(G))]
        report('ROW ' + ''.join(map(str, a_)) + f' ({len(G)} states)', fromF(Vr), homog)
print(f'time {time.time() - t0:.1f}s')
print('PASS')
