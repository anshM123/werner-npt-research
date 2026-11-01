"""Exact linear data at n = 4 for the lifted (all-ratio) problem: the swap-test functionals fu, fv, fuv, the objectives
    c+ = (fu + fv + 2 fuv) / 16   (q_4(P_u + P_v), the exact4 objective),
    c- = (fu + fv - 2 fuv) / 16   (q_4(P_u - P_v), the exact4m objective),
    d  = (fu - fv) / 16           (q_4(P_u) - q_4(P_v)),
so that q_4(a P_u + b P_v) = (a^2 fu + b^2 fv + 2 a b fuv) . x / 16 = (1/4) [(a+b)^2 c+ + (a-b)^2 c- + 2 (a^2 - b^2) d] . x,
and the equalities E x = e (120 moment functionals + ORTH/MARG identity rows, exact RREF), exactly as in
../exact4/build_lin4.py (same code path: fu = Tr[X_uuuu W01], fv = Tr[X_vvvv W01], fuv = Tr[X_uuvv W12],
W_ij = (x)_s (2 F_ij - I)).
Checks (exact): c+ equals the certified objective of ../exact4/lin4.pkl, c- equals the certified objective of
../exact4m/lin4m.pkl (and its cplus equals c+), E and e equal the saved ones in both files.
Saves linL.pkl.   usage: python build_linL.py"""
import sys, time, pickle, itertools, math, numpy as np, flint
from gmpy2 import mpq
from core4 import *
from engine4 import *

P_SEL = 1048573          # < 2^20: float64 products of residues over <= 4096 terms stay below 2^53

def independent_rows_modp(M, p=P_SEL, chunk=128):
    """indices of a maximal set of rows of the integer matrix M (int64) independent modulo p"""
    Mp = np.mod(M, p).astype(np.float64)
    nr, nc = Mp.shape
    B = np.zeros((0, nc)); piv = []; chosen = []
    for c0 in range(0, nr, chunk):
        C = Mp[c0:c0 + chunk].copy()
        if piv:
            C = np.mod(C - np.mod(C[:, piv] @ B, p), p)
        for r in range(C.shape[0]):
            row = C[r]
            nz = np.nonzero(row)[0]
            if nz.size == 0: continue
            j = int(nz[0]); inv_ = pow(int(row[j]), -1, p)
            row = np.mod(row * inv_, p)
            if B.shape[0]:
                B = np.mod(B - np.outer(B[:, j], row), p)
            C[r + 1:] = np.mod(C[r + 1:] - np.outer(C[r + 1:, j], row), p)
            B = np.vstack([B, row]); piv.append(j); chosen.append(c0 + r)
    return chosen

def main():
    t0 = time.time(); n = 4
    A, B1, A3 = site_algebras()
    vs = VarSpace(A, n); N = vs.N
    def functional_row(c, w):
        Ti, d = functional_map_int(A, w)
        row = [mpq(0)] * N
        for key in vs.order:
            if key[0] != c: continue
            g = vs.groups[key]
            level = tree(vs.pair_float(key), key[1], Ti, 1)
            acc = assemble(level, g['orbit'], (0,) * n)
            if acc is None: continue
            check_exact(acc, g['Xmax'] * rowsum_bound(Ti) ** n * len(g['orbit']))
            vals = acc.reshape(g['K'])
            for k, i in enumerate(g['idx']): row[i] = mpq(int(vals[k]) * g['mlam'], d ** n * g['dX'])
        return row
    w01 = PERM[(1, 0, 2, 3)] * 2 - np.eye(81, dtype=np.int64); w12 = PERM[(0, 2, 1, 3)] * 2 - np.eye(81, dtype=np.int64)
    fu = functional_row('uuuu', w01); fv = functional_row('vvvv', w01); fuv = functional_row('uuvv', w12)
    s = mpq(1, 2 ** n)
    cplus = [(fu[i] + fv[i] + 2 * fuv[i]) * s for i in range(N)]
    cminus = [(fu[i] + fv[i] - 2 * fuv[i]) * s for i in range(N)]
    d = [(fu[i] - fv[i]) * s for i in range(N)]
    mrows = []; mrhs = []
    for cf in FAMS:
        for tau in S4:
            val = 1; seen = set()
            for i in range(4):
                if i in seen: continue
                j = i; letters = set()
                while j not in seen: seen.add(j); letters.add(cf[j]); j = tau[j]
                if len(letters) > 1: val = 0
            mrows.append(functional_row(cf, PERM[tau])); mrhs.append(mpq(val))
    print(f"[lin] fu, fv, fuv and {len(mrows)} moment rows; {time.time() - t0:.1f}s", flush=True)
    blocks = []; names = []
    for name, terms in identity_defs():
        maps = [(coef, cf) + site_map_int(A, A3, phi) for coef, cf, phi in terms]
        Lrow = 1
        for _, _, Ti, dT in maps: Lrow = Lrow * dT ** n // math.gcd(Lrow, dT ** n)
        mus = sorted_patterns(len(A3.blocks), n)
        Rs = {mu: np.zeros((int(np.prod([A3.blocks[q]['f'] for q in mu])) ** 2, N), dtype=np.int64) for mu in mus}
        for coef, cf, Ti, dT in maps:
            for key in vs.order:
                if key[0] != cf: continue
                g = vs.groups[key]
                level = tree(vs.pair_float(key), key[1], Ti, len(A3.blocks))
                bnd = g['Xmax'] * rowsum_bound(Ti) ** n * len(g['orbit'])
                assert bnd * (Lrow // dT ** n) < 2 ** 62
                for mu in mus:
                    acc = assemble(level, g['orbit'], mu)
                    if acc is None: continue
                    check_exact(acc, bnd)
                    gs = tuple(A3.blocks[q]['f'] for q in mu); Dm = int(np.prod(gs))
                    Y = pair_to_mat(acc, gs).reshape(g['K'], Dm * Dm).astype(np.int64)
                    Rs[mu][:, g['idx']] += coef * (Lrow // dT ** n) * Y.T
        for mu in mus:
            if np.any(Rs[mu]):
                blocks.append(Rs[mu]); names.append((name, mu))
    Rt = np.vstack(blocks)
    sel = independent_rows_modp(Rt)
    print(f"[lin] identity rows: {Rt.shape[0]} raw rows, rank mod {P_SEL}: {len(sel)}; {time.time() - t0:.1f}s", flush=True)
    dXcol = [vs.groups[vs.order[vs.gid[i]]]['dX'] for i in range(N)]
    nrow = len(mrows) + len(sel)
    ent = []
    for r, rhs in zip(mrows, mrhs):
        for v in r: ent.append(fq(int(v.numerator), int(v.denominator)))
        ent.append(fq(int(rhs.numerator), int(rhs.denominator)))
    for k in sel:
        row = Rt[k]
        for i in range(N): ent.append(fq(int(row[i]), dXcol[i]))
        ent.append(fq(0))
    Mq = flint.fmpq_mat(nrow, N + 1, ent)
    RR, rk = Mq.rref()
    E = [[mpq(int(RR[i, j].p), int(RR[i, j].q)) for j in range(N)] for i in range(rk)]
    e = [mpq(int(RR[i, N].p), int(RR[i, N].q)) for i in range(rk)]
    for i in range(rk):
        assert any(v != 0 for v in E[i]), "inconsistent equalities"
    Iq = flint.fmpq_mat(len(sel), N, [fq(int(Rt[k][i]), dXcol[i]) for k in sel for i in range(N)])
    rkI = Iq.rank()
    nbits = max(max(int(abs(v.numerator)).bit_length(), int(v.denominator).bit_length()) for row in E for v in row)
    print(f"[lin] equalities: {rk} independent (identity rows alone: exact rank {rkI}); max E entry {nbits} bits; "
          f"{time.time() - t0:.1f}s", flush=True)
    # ---------------------------------------------------------------- exact comparison with the certified data
    l4 = pickle.load(open('../exact4/lin4.pkl', 'rb')); l4m = pickle.load(open('../exact4m/lin4m.pkl', 'rb'))
    ok_cp = all(cplus[i] == l4['c'][i] for i in range(N)) and len(l4['c']) == N
    ok_cm = all(cminus[i] == l4m['c'][i] for i in range(N)) and l4m['t'] == -1
    ok_cpm = all(cplus[i] == l4m['cplus'][i] for i in range(N))
    ok_E = (len(E) == len(l4['E']) == len(l4m['E']) and
            all(E[j][i] == l4['E'][j][i] == l4m['E'][j][i] for j in range(len(E)) for i in range(N)))
    ok_e = all(e[j] == l4['e'][j] == l4m['e'][j] for j in range(len(e)))
    ok_m = all(mrows[j][i] == l4['mrows'][j][i] for j in range(len(mrows)) for i in range(N))
    ok_rel = all(cplus[i] - cminus[i] == fuv[i] * mpq(4, 2 ** n) for i in range(N))
    nd = sum(1 for v in d if v != 0)
    print(f"[lin] exact comparison: c+ == ../exact4/lin4.pkl c: {ok_cp}; c- == ../exact4m/lin4m.pkl c: {ok_cm}; "
          f"c+ == lin4m cplus: {ok_cpm}; E identical in both: {ok_E}; e identical: {ok_e}; moment rows identical: {ok_m}; "
          f"c+ - c- == fuv/4: {ok_rel}; d has {nd} nonzero coordinates", flush=True)
    assert ok_cp and ok_cm and ok_cpm and ok_E and ok_e and ok_m and ok_rel
    with open('linL.pkl', 'wb') as fh:
        pickle.dump(dict(N=N, fu=fu, fv=fv, fuv=fuv, cplus=cplus, cminus=cminus, d=d, E=E, e=e, mrows=mrows, mrhs=mrhs,
                         ident_rank=rkI, n_ident_raw=Rt.shape[0]), fh)
    print(f"saved linL.pkl; {time.time() - t0:.1f}s", flush=True)

if __name__ == "__main__":
    main()
