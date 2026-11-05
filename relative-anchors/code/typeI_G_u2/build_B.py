"""Condition (B_2) as polynomials: lambda_min(Ghat_tt) > H (1 - s_t/sigma) for u = 2.

Normalisation sigma + S = 1 (sigma = s_T, S = s_U): s_t = th (1 - mu), sigma = 1 - mu, s_3 = mu la, s_4 = mu (1 - la),
(th, mu, la) in [0,1]^3.  With G_z = G00 in the coordinates zeta_{t,h} = w_t w_h z_{t,h} (closed_u2.pkl) and
Ghat = D G_zeta D (D^2 = I + w_U w_U^T/sigma on the zeta part), Ghat - h I > 0 iff
    N = G_z - h diag(1, s_t (diag(s_U) - s_U s_U^T/(sigma + S))) > 0,     h = H (1 - th).
Multiplying by the common denominator den > 0 gives Nt = den N with polynomial entries; Nt > 0 iff its three
leading principal minors are positive.  Each minor is divided by the largest powers of th, mu, la, 1 - la that divide it
(positive on the open box).  Output: build_B.pkl = {H: [(name, dict, degs, removed factors)]}."""
import pickle, sys, flint, time
from fractions import Fraction as Fr

ctx5 = flint.fmpq_mpoly_ctx.get(('s0', 's1', 's2', 's3', 's4'), 'lex')
ctx = flint.fmpq_mpoly_ctx.get(('th', 'mu', 'la'), 'lex')
th, mu, la = ctx.gens()
one = ctx.from_dict({(0, 0, 0): 1})


def load():
    D = pickle.load(open('closed_u2.pkl', 'rb'))
    sig = 1 - mu
    st = th * sig
    sub = (st, (sig - st) * flint.fmpq(1, 2), (sig - st) * flint.fmpq(1, 2), mu * la, mu * (1 - la))   # s1, s2 enter via s1 + s2 only
    num = [[ctx5.from_dict(D['G00'][i][j][0]).compose(*sub, ctx=ctx) for j in range(3)] for i in range(3)]
    den = [[ctx5.from_dict(D['G00'][i][j][1]).compose(*sub, ctx=ctx) for j in range(3)] for i in range(3)]
    d0 = den[0][0]
    assert all(den[i][j] == d0 for i in range(3) for j in range(3))
    return num, d0


def strip(p):
    """divide p by the largest powers of th, mu, la, (1 - la) dividing it."""
    removed = []
    for nm, f in (('th', th), ('mu', mu), ('la', la), ('1-la', 1 - la), ('1-th', 1 - th), ('1-mu', 1 - mu)):
        e = 0
        while True:
            q, r = divmod(p, f)
            if r != 0:
                break
            p = q; e += 1
        if e:
            removed.append((nm, e))
    return p, removed


def build(H):
    num, den = load()
    sig = 1 - mu; st = th * sig
    sU = [mu * la, mu * (1 - la)]
    h = flint.fmpq(H.numerator, H.denominator) * (1 - th)
    Oinv = [[one, 0 * one, 0 * one], [0 * one, st * (sU[0] - sU[0] * sU[0]), -st * sU[0] * sU[1]],
            [0 * one, -st * sU[0] * sU[1], st * (sU[1] - sU[1] * sU[1])]]   # sigma + S = 1
    Nt = [[num[i][j] - h * den * Oinv[i][j] for j in range(3)] for i in range(3)]
    m1 = Nt[0][0]
    m2 = Nt[0][0] * Nt[1][1] - Nt[0][1] * Nt[1][0]
    m3 = (Nt[0][0] * (Nt[1][1] * Nt[2][2] - Nt[1][2] * Nt[2][1]) - Nt[0][1] * (Nt[1][0] * Nt[2][2] - Nt[1][2] * Nt[2][0])
          + Nt[0][2] * (Nt[1][0] * Nt[2][1] - Nt[1][1] * Nt[2][0]))
    out = []
    for name, p in (('m1', m1), ('m2', m2), ('m3', m3)):
        q, rem = strip(p)
        out.append((name, q, rem))
    return out, den


if __name__ == '__main__':
    t0 = time.time()
    res = {}
    for H in (Fr(2, 5), Fr(9, 20)):
        out, den = build(H)
        lst = []
        for name, q, rem in out:
            dct = {m: Fr(int(c.p), int(c.q)) for m, c in q.to_dict().items()}
            lst.append((name, dct, q.degrees(), rem))
            print(f'H = {H}: {name}: degrees (th, mu, la) = {q.degrees()}, terms {len(q)}, removed {rem}', flush=True)
        res[H] = lst
    dd = {m: Fr(int(c.p), int(c.q)) for m, c in den.to_dict().items()}
    res['den'] = dd
    pickle.dump(res, open('build_B.pkl', 'wb'))
    print(f'time {time.time() - t0:.1f} s')
    print('PASS')
