"""Exact tools: load gauge-fixed blocks as polynomial matrices Nt = den * M_g, monomial coefficient matrices, an exact
PSD test with a rational witness, and helpers for coefficient signs."""
import pickle
from fractions import Fraction as Fr
import flint


def load_blocks(k):
    """returns list of dicts with keys: name, k, j, kind, T, labels (kept), states (kept), names (variables),
    Mg (list of lists of (num fmpq_mpoly, den fmpq_mpoly)), ctx."""
    out = []
    if k in (2, 3):
        D = pickle.load(open(f'blocks_k{k}.pkl', 'rb'))
        names = tuple(D['names'])
        ctx = flint.fmpq_mpoly_ctx.get(names, 'lex')
        for b in D['blocks']:
            keep = b['keep']
            lab = b['label']
            T = tuple(t for t, x in enumerate(lab[1:]) if x == 1)
            Mg = [[(ctx.from_dict(b['M'][a][c][0]), ctx.from_dict(b['M'][a][c][1])) for c in keep] for a in keep]
            out.append(dict(name=f'k{k}_j{b["j"]}_{lab}', k=k, j=b['j'], kind=b['kind'], label=lab, T=T,
                            labels=[b['labels'][a] for a in keep], states=[b['states'][a] for a in keep],
                            names=names, Mg=Mg, ctx=ctx))
    elif k == 4:
        D = pickle.load(open('blocks_sym_k4.pkl', 'rb'))
        names = tuple(D['names'])
        ctx = flint.fmpq_mpoly_ctx.get(names, 'lex')
        for (j, T), b in D['blocks'].items():
            keep = b['keep']
            kind = 'L1' if len(T) == 0 else 'I'
            Mg = [[(ctx.from_dict(b['M'][a][c][0]), ctx.from_dict(b['M'][a][c][1])) for c in keep] for a in keep]
            out.append(dict(name=f'k4_j{j}_T{"".join(map(str, T)) or "-"}', k=4, j=j, kind=kind, label=None, T=T,
                            labels=[b['labels'][a] for a in keep], states=[b['states'][a] for a in keep],
                            names=names, Mg=Mg, ctx=ctx))
    return out


def lcm_poly(ps):
    den = None
    for d in ps:
        if den is None:
            den = d
        else:
            g = den.gcd(d)
            den = den * (d / g)
    return den


def poly_matrix(blk, mult=None):
    """Nt = den * M_g (* mult), den = lcm of entry denominators normalised to positive leading coefficient."""
    Mg = blk['Mg']
    n = len(Mg)
    den = lcm_poly([Mg[a][c][1] for a in range(n) for c in range(n)])
    if den.leading_coefficient() < 0:
        den = -den
    Nt = [[Mg[a][c][0] * (den / Mg[a][c][1]) for c in range(n)] for a in range(n)]
    if mult is not None:
        Nt = [[x * mult for x in row] for row in Nt]
    return den, Nt


def coeff_mats(Nt):
    """dict monomial -> dense Fraction matrix of coefficients."""
    n = len(Nt)
    mons = {}
    for a in range(n):
        for c in range(n):
            for m, v in Nt[a][c].to_dict().items():
                mons.setdefault(m, {})[(a, c)] = Fr(int(v.p), int(v.q))
    out = {}
    for m, ent in mons.items():
        C = [[Fr(0)] * n for _ in range(n)]
        for (a, c), v in ent.items():
            C[a][c] = v
        out[m] = C
    return out


def psd_exact(C):
    """(True, None) if the symmetric rational matrix C is PSD, else (False, x) with x^T C x < 0 (x rational)."""
    n = len(C)
    A = [row[:] for row in C]
    active = list(range(n))
    piv = []  # (i, row dict over the then-active set)
    while True:
        cand = [i for i in active if A[i][i] > 0]
        if not cand:
            break
        i = cand[0]
        rest = [r for r in active if r != i]
        row = {r: A[i][r] for r in rest}
        piv.append((i, A[i][i], row))
        for r in rest:
            if A[r][i] != 0:
                f = A[r][i] / A[i][i]
                for l in rest:
                    if A[i][l] != 0:
                        A[r][l] -= f * A[i][l]
        active = rest
    z = None
    for l in active:
        if A[l][l] < 0:
            z = {l: Fr(1)}
            break
    if z is None:
        for a_ in active:
            for b_ in active:
                if a_ < b_ and A[a_][b_] != 0:
                    z = {a_: Fr(1), b_: Fr(-1) if A[a_][b_] > 0 else Fr(1)}
                    break
            if z is not None:
                break
    if z is None:
        return True, None
    x = {r: Fr(0) for r in range(n)}
    for r, v in z.items():
        x[r] = v
    for (i, d, row) in reversed(piv):
        x[i] = -sum(row[r] * x[r] for r in row) / d
    xv = [x[r] for r in range(n)]
    val = sum(xv[a] * C[a][b] * xv[b] for a in range(n) for b in range(n))
    assert val < 0, val
    return False, xv


def sign_profile(p):
    """'+', '-', 'mixed' or '0' for a polynomial."""
    if p.is_zero():
        return '0'
    cs = p.coeffs()
    if all(c > 0 for c in cs):
        return '+'
    if all(c < 0 for c in cs):
        return '-'
    return 'mixed'
