"""Condition (A_2) as a polynomial: ||Ghat_tt'||_F <= H sqrt(s_t s_t') / sigma for u = 2.

Normalisation sigma + S = 1: s_t = th (1 - mu), s_t' = ph (1 - mu), s_v = (1 - th - ph)(1 - mu) (the other T-copies),
s_3 = mu la, s_4 = mu (1 - la).  With G01 = Nm/den (z-coordinates, closed_u2.pkl) and Omega_t = Z^-1 D^2 Z^-1 =
diag(1, (diag(1/s_U) + 1 1^T/sigma)/s_t), ||Ghat_tt'||_F^2 = tr(Omega_t G01 Omega_t' G01^T).  Put
Om~_t = s_t s_3 s_4 sigma Omega_t = diag(s_t s_3 s_4 sigma, sigma diag(s_4, s_3) + s_3 s_4 1 1^T) (polynomial).  Then
(A_2) <=> FA := H^2 s_t^2 s_t'^2 (s_3 s_4)^2 den^2 - tr(Om~_t Nm Om~_t' Nm^T) >= 0."""
import pickle, sys, flint, time
from fractions import Fraction as Fr

ctx5 = flint.fmpq_mpoly_ctx.get(('s0', 's1', 's2', 's3', 's4'), 'lex')
ctx = flint.fmpq_mpoly_ctx.get(('th', 'ph', 'mu', 'la'), 'lex')
th, ph, mu, la = ctx.gens()


def strip(p):
    removed = []
    for nm, f in (('th', th), ('ph', ph), ('mu', mu), ('la', la), ('1-la', 1 - la), ('1-mu', 1 - mu), ('1-th-ph', 1 - th - ph)):
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
    D = pickle.load(open('closed_u2.pkl', 'rb'))
    sig = 1 - mu
    s0, s1 = th * sig, ph * sig
    s3, s4 = mu * la, mu * (1 - la)
    sub = (s0, s1, sig - s0 - s1, s3, s4)
    Nm = [[ctx5.from_dict(D['G01'][i][j][0]).compose(*sub, ctx=ctx) for j in range(3)] for i in range(3)]
    dn = [[ctx5.from_dict(D['G01'][i][j][1]).compose(*sub, ctx=ctx) for j in range(3)] for i in range(3)]
    den = dn[0][0]
    assert all(dn[i][j] == den for i in range(3) for j in range(3))
    zero = 0 * th
    def Om(st):
        return [[st * s3 * s4 * sig, zero, zero], [zero, sig * s4 + s3 * s4, s3 * s4], [zero, s3 * s4, sig * s3 + s3 * s4]]
    O1, O2 = Om(s0), Om(s1)
    # tr(O1 Nm O2 Nm^T) = sum_{i,j,k,l} O1[i][j] Nm[j][k] O2[k][l] Nm[i][l]
    A = [[sum((Nm[j][k] * O2[k][l] for k in range(3)), zero) for l in range(3)] for j in range(3)]   # Nm O2
    tr = zero
    for i in range(3):
        for j in range(3):
            if O1[i][j] == 0:
                continue
            for l in range(3):
                tr += O1[i][j] * A[j][l] * Nm[i][l]
    Hq = flint.fmpq(H.numerator, H.denominator)
    FA = Hq * Hq * s0 * s0 * s1 * s1 * s3 * s3 * s4 * s4 * den * den - tr
    return FA, den


if __name__ == '__main__':
    t0 = time.time()
    H = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2, 5)
    FA, den = build(H)
    q, rem = strip(FA)
    print(f'H = {H}: FA degrees (th, ph, mu, la) = {q.degrees()}, terms {len(q)}, removed {rem} ({time.time() - t0:.1f} s)', flush=True)
    dct = {m: Fr(int(c.p), int(c.q)) for m, c in q.to_dict().items()}
    dd = {m: Fr(int(c.p), int(c.q)) for m, c in den.to_dict().items()}
    pickle.dump(dict(H=H, FA=dct, degs=q.degrees(), removed=rem, den=dd), open('build_A.pkl', 'wb'))
    print('PASS')
