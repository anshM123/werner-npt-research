"""Check of Section 10 (exact): N~/Den_c against det A computed directly at points far outside the grid x1 = 0..B.

At each point (u1, u2, sigma, s1) the closed form tl_inv_gen_T11/T22 is evaluated exactly (rational numbers, no
polynomial algebra), G_y = N^T G_z N, kappa and A = G_y - H (r/sigma) diag(kappa) are formed as in check_fN.py,
part 2 (an implementation independent of em3lib.py), and det A and the three pivots are computed over Q.  This is
compared with N~(x1, x2, r, s)/Den_c(x1, x2, r, s) at x_i = u_i - 2, r = sigma - s_t, from the interpolated
polynomial N~ (em3_out/N_KA.pkl) and the factors of Den_c (em3_out/setup_KA.pkl).  Every point has x1 outside the
grid {0, ..., B} (x1 = 77 is the first value above it); one point has u1 and u2 non-integer (the identity
det A = N~/Den_c is one of rational functions).
python em3points.py"""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint

Q = flint.fmpq
H = Q(1, 2)
OUT = 'em3_out'
CTX4 = flint.fmpq_mpoly_ctx.get(('x1', 'x2', 'r', 's'), 'lex')
POINTS = [(Q(300), Q(17), Q(2, 7), Q(5, 3)),
          (Q(1000), Q(1000), Q(1000), Q(1, 1000)),
          (Q(79), Q(5000), Q(1, 1000), Q(50)),
          (Q(2003, 7), Q(41, 3), Q(7, 5), Q(3)),
          (Q(100000), Q(23), Q(3, 11), Q(1, 7))]          # (u1, u2, s, r)


def load_closed(md):
    D = pickle.load(open('tl_inv_gen_%s.pkl' % md, 'rb'))
    ctx = flint.fmpq_mpoly_ctx.get(tuple(D['names']), 'lex')
    P = lambda d: ctx.from_dict({m: Q(p, q) for m, (p, q) in d.items()})
    return D['names'], {k: (P(n), P(d)) for k, (n, d) in D['G'].items()}


def direct(names, G, u1, u2, s1, st, rr):
    """det A and the three pivots at one point, from the closed form evaluated over Q"""
    sig = st + rr
    vals = {'u1': u1, 'u2': u2, 'sig': sig, 's1': s1}
    args = [vals[x] for x in names]
    Gz = [[G[(min(i, j), max(i, j))][0](*args) / G[(min(i, j), max(i, j))][1](*args) for j in range(3)]
          for i in range(3)]
    S1, S2 = u1 * s1, u2
    S = S1 + S2
    Nm = [[Q(1), Q(0), Q(0)], [Q(0), 1 / S, S2 / S], [Q(0), 1 / S, -S1 / S]]
    Gy = [[sum((Nm[k][i] * Gz[k][l] * Nm[l][j] for k in range(3) for l in range(3)), Q(0)) for j in range(3)]
          for i in range(3)]
    kap = [Q(1), st * sig / (S * (sig + S)), st * S1 * S2 / S]
    M = [[Gy[i][j] - (H * (rr / sig) * kap[i] if i == j else 0) for j in range(3)] for i in range(3)]
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    p1 = M[0][0]
    p2 = M[1][1] - M[0][1] * M[1][0] / M[0][0]
    return det, (p1, p2, det / (p1 * p2))


def main():
    t0 = time.time()
    ok = True
    for ka, md in (('s', 'T11'), ('1', 'T22')):
        st_ = pickle.load(open(os.path.join(OUT, 'setup_%s.pkl' % ka), 'rb'))
        Nd = pickle.load(open(os.path.join(OUT, 'N_%s.pkl' % ka), 'rb'))
        Npoly = CTX4.from_dict({k: Q(p, q) for k, (p, q) in Nd['N'].items()})
        facs = [CTX4.from_dict({m: Q(p, q) for m, (p, q) in fd.items()}) for fd in st_['fac_dicts']]
        names, G = load_closed(md)
        print('s_t = %s: N~ with %d terms loaded [%.0fs]' % (ka, len(Npoly), time.time() - t0), flush=True)
        for (u1, u2, s, rr) in POINTS:
            st = s if ka == 's' else Q(1)
            pt = (u1 - 2, u2 - 2, rr, s)
            dc = Q(1)
            for f, e in zip(facs, st_['dexp']):
                dc *= f(*pt) ** e
            val = Npoly(*pt) / dc
            det, piv = direct(names, G, u1, u2, s, st, rr)
            good = (val == det) and all(p > 0 for p in piv)
            ok &= good
            print('  (u1, u2, s, r) = (%s, %s, %s, %s): N~/Den_c = det A: %s; det A = %.6e; pivots %s'
                  % (u1, u2, s, rr, val == det, float(det), ', '.join('%.4e' % float(p) for p in piv)), flush=True)
    print('all equal, all pivots positive: %s [%.0fs]' % (ok, time.time() - t0))
    return ok


if __name__ == '__main__':
    print('PASS' if main() else 'FAIL')
