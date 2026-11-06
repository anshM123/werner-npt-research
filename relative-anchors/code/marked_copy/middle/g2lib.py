"""The generic two-row family [q0 >= 4, t >= 4, i >= 7] of the letter-1 block with one marked copy, by the
determinant route: the gauge-fixed form Q (10 x 10, U10 and U11 deleted) over Q(q0, t, i, s), the multiple D* of the
denominator of det Q from the factored entry denominators, and rigorous degree bounds for T = D* det Q.

D*: for every irreducible factor f of an entry denominator, with m_f(a, b) its multiplicity in den Q_ab, the
exponent of f in D* is max over permutations p of sum_a m_f(a, p(a)).  Then D* * prod_a Q_{a,p(a)} is a polynomial
for every p, so T = D* det Q is a polynomial.
Degree bound in a variable v: deg_v T <= deg_v D* + max_p sum_a (deg_v num Q_{a,p(a)} - deg_v den Q_{a,p(a)}) over
the permutations with all Q_{a,p(a)} nonzero (each term D* prod Q_{a,p(a)} is a polynomial of at most that degree).
The maxima over permutations are computed exactly by dynamic programming over subsets."""
import time
import flint
import famlib
famlib.SMIN['q0'] = 4
famlib.SMIN['t'] = 4
famlib.SMIN['i'] = 7
import splitlib
from prove_fam import one_signed, fam_kernel, det_small
from prove_split import sanity

GAUGE = ('U10', 'U11')
VARS = ('q0', 't', 'i', 's')


def build(check=True, log=print):
    """the 4-variable family form; with check=True the steps of prove_gen.prove_gen before the elimination (sanity
    against mlib, entry denominators one-signed, split block-diagonal, kernel vectors annihilated, gauge determinant
    one-signed).  Returns (C, keep, M) with M[a, b] the entries of the gauge-fixed two-row form."""
    t0 = time.time()
    C = famlib.FamCalc(None, None, None)
    names, Q = famlib.reduced(C)
    log(f'  form computed ({time.time()-t0:.1f}s), {len(names)} vectors')
    if check:
        bad = sanity(C, names, Q, None, None, None, 1)
        assert not bad, bad
        log(f'  sanity against mlib done ({time.time()-t0:.1f}s)')
        n = len(names)
        for a in range(n):
            for b in range(n):
                if Q[a][b] != 0:
                    assert one_signed(Q[a][b].d, C) != 0, ('entry denominator', names[a], names[b])
    vecs, split = splitlib.basis(names, C.q0, C.m, C.i, C.one)
    d = C.m + 1 + C.i
    assert one_signed(d.n, C) == 1 and one_signed(d.d, C) == 1
    two, three, M = splitlib.congruence(names, Q, vecs, C.zero)
    for a in two:
        for b in three:
            assert M[a, b] == 0, ('split', a, b)
    if check:
        fk = fam_kernel(names, C)
        kv = []
        for lab, v, dl in fk:
            w = splitlib.to_new(v, split)
            y0 = lab[-1]
            assert w is not None and w.get('U1' + y0) == 1, lab
            assert all(nm in two for nm in w), lab
            kv.append((lab, w, 'U1' + y0))
        for lab, v, _ in kv:
            for a in two:
                s = C.zero
                for nm, cf in v.items():
                    s = s + M[a, nm] * cf
                assert s == 0, ('kappa', lab, a)
        gauge = [dl for _, _, dl in kv]
        assert tuple(gauge) == GAUGE, gauge
        sub = [[v.get(c, C.zero) for c in gauge] for _, v, _ in kv]
        det = det_small(sub, C)
        assert det != 0 and one_signed(det.n, C) != 0 and one_signed(det.d, C) != 0
        log(f'  split, kernel vectors and gauge checked ({time.time()-t0:.1f}s)')
    keep = [nm for nm in two if nm not in GAUGE]
    return C, keep, M


def maxperm(W):
    """max over permutations p of sum_a W[a][p(a)] (W[a][b] = None: forbidden); exact DP over subsets"""
    n = len(W)
    NEG = None
    best = {0: 0}
    for a in range(n):
        nb = {}
        for mask, val in best.items():
            for b in range(n):
                if mask >> b & 1 or W[a][b] is None:
                    continue
                m2 = mask | (1 << b)
                v = val + W[a][b]
                if m2 not in nb or v > nb[m2]:
                    nb[m2] = v
        best = nb
    return best.get((1 << n) - 1, NEG)


def dstar(C, keep, M):
    """the irreducible factors of the entry denominators and the exponents of D*; returns (Dstar, [(f, e)])"""
    facs = {}
    mult = {}
    for a in keep:
        for b in keep:
            e = M[a, b]
            if e == 0:
                continue
            c, fl = e.d.factor()
            for f, m in fl:
                key = str(f)
                facs[key] = f
                mult[(key, a, b)] = m
    out = []
    D = C.K.ctx.from_dict({(0,) * 4: 1})
    for key, f in facs.items():
        W = [[(mult.get((key, a, b), 0) if M[a, b] != 0 else None) for b in keep] for a in keep]
        e = maxperm(W)
        out.append((f, e))
        D = D * f ** e
    return D, out


def degree_bounds(C, keep, M, D):
    """rigorous upper bounds for deg_v (D* det Q), v = q0, t, i, s"""
    out = {}
    for iv, v in enumerate(VARS):
        W = [[(M[a, b].n.degrees()[iv] - M[a, b].d.degrees()[iv]) if M[a, b] != 0 else None for b in keep]
             for a in keep]
        out[v] = D.degrees()[iv] + maxperm(W)
    return out
