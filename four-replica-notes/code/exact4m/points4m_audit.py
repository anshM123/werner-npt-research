"""Audit copy of points4m.py: adds kinds near_<base> (a zero scaled by 4 plus a small integer perturbation, then v
orthogonalised against u exactly).  Everything else is unchanged.

Exact points of the relaxation at n = 4, t = -1: zeros of q_4(P_u - P_v) and entangled test points.

  vectors   integer u, v in Z^81 with <u, v> = 0 and |u| = |v| (zeros) or u = p, v = (p.p) q - (p.q) p (test points);
            zero families (all real, integer):
              k2     u = a (x) b (x) psi, v = a' (x) b' (x) psi; (a, a'), (b, b') two rows of integer Cayley orthogonal
                     matrices on sites 1, 2; psi in Z^9 random with full Schmidt rank on sites 3, 4;
              k4     u = a (x) b (x) c (x) d, v = a' (x) b' (x) c' (x) d' (orthogonal on all four sites);
              inner  u = a (x) phi, v = a' (x) chi, (phi, chi) = (r1 p + r2 q, -r2 p + r1 q), (r1, r2) a Pythagorean pair,
                     p = b (x) c (x) d, q = b' (x) c' (x) d'  (q_4(P_{a phi} - P_{a' chi}) = q_3(P_phi + P_chi) / 2 and
                     P_phi + P_chi is proportional to the projector onto span(p, q), a t = +1 zero at three copies);
              sf     spin flip: u supported on the qubit sector {0,1}^4, u_x in {-1, 0, 1}, u_0000 = 1, u_1111 = -S with
                     S = sum over the other 7 pairs {x, xbar} of (-1)^w(x) u_x u_xbar, so that u^T Y u = 0 for
                     Y|x> = (-1)^w(x) |xbar> (= sigma_y^(x)4, real); v = Y u (then |u| = |v|, <u, v> = 0);
  moments   m[sigma] = Tr[(x)_r |c_r><c_r| (sigma_1 (x) ... (x) sigma_4)] (integer, unnormalised) modulo 27-bit
            primes (pairwise einsum with an asserted bound: summed terms * p^2 < 2^63) and CRT with the bound
            prod_r ||c_r||_1^2 (as in ../relaxations/exact_coreN.moments_exact);
  hats      exact integer Fourier transform per site (F[(l,i,j), sigma] = f_l/(24 m_l) rho_l(sigma^-1)[i,j]),
            site symmetrisation, exact X_hat per sorted pattern;
  x         exact projection onto the parameters with exactness checks (VarSpace.x_from_hats).
q_4 is computed directly from the 81 x 81 matrix P = P_u - P_v (and, as a consistency check of the objective pieces,
from P_u + P_v against the t = +1 objective cplus).
Saves pts/<name>.pkl with x (list of (num, den)), the vectors, and the exact checks.
usage: python points4m.py <kind> <seed> <name>     kind in k2, k4, inner, sf, ent"""
import sys, os, time, pickle, itertools, math, numpy as np, flint
from gmpy2 import mpq
from core4 import *
from engine4 import *
import exact_coreN as XN

import gmpy2
def _primes_below(B, k):
    out = []; q = B
    while len(out) < k:
        q -= 1
        if gmpy2.is_prime(q): out.append(q)
    return out
PRIMES = _primes_below(2 ** 27, 16)      # 243 * (2^27)^2 < 2^63: the largest pairwise contraction at n = 4 sums 3^5 terms

def _contract_mod_safe(ops, subs, out, p):
    ops = [o % p for o in ops]; subs = list(subs)
    while len(ops) > 1:
        best = None
        for i in range(len(ops)):
            for j in range(i + 1, len(ops)):
                sh = set(subs[i]) & set(subs[j])
                if sh:
                    others = set(''.join(subs[k] for k in range(len(ops)) if k not in (i, j))) | set(out)
                    keep = ''.join(ch for ch in subs[i] + subs[j] if (ch not in sh or ch in others))
                    keep = ''.join(dict.fromkeys(keep))
                    size = 3 ** len(keep)
                    if best is None or size < best[0]: best = (size, i, j, keep)
        if best is None:
            i, j = 0, 1; keep = ''.join(dict.fromkeys(subs[0] + subs[1]))
        else:
            _, i, j, keep = best
        summed = set(subs[i] + subs[j]) - set(keep)
        assert (3 ** len(summed)) * (p - 1) ** 2 < 2 ** 63, "int64 overflow risk in contraction"
        res = np.einsum(subs[i] + ',' + subs[j] + '->' + keep, ops[i], ops[j]) % p
        ops = [ops[k] for k in range(len(ops)) if k not in (i, j)] + [res]
        subs = [subs[k] for k in range(len(subs)) if k not in (i, j)] + [keep]
    res = ops[0]
    if subs[0] != out:
        assert set(subs[0]) == set(out)
        res = np.einsum(subs[0] + '->' + out, res) % p
    return res

def moments_modp(Ls, n, p):
    letters = [chr(c) for c in list(range(65, 91)) + list(range(97, 123))]
    a = letters[:4]; b = letters[4:8]
    X = [[letters[8 + r * n + s] for s in range(n)] for r in range(4)]
    T = [np.array([int(v) % p for v in L], dtype=np.int64).reshape((3,) * n) for L in Ls]
    PTf = np.array([PERM[s].T.ravel() for s in S4], dtype=np.int64)
    out = np.zeros((24,) * n, dtype=np.int64)
    for head in itertools.product(range(24), repeat=n - 1):
        sig = [S4[k] for k in head]
        subs = []; ops = []
        for r in range(4):
            subs.append(''.join(X[r][s] for s in range(n - 1)) + a[r]); ops.append(T[r])
            subs.append(''.join(X[inv(sig[s])[r]][s] for s in range(n - 1)) + b[r]); ops.append(T[r])
        M = _contract_mod_safe(ops, subs, ''.join(a) + ''.join(b), p).reshape(-1)
        out[head] = (PTf @ M) % p          # each row of PTf has 81 ones: sums < 81 p
    return out

def moments_int(Ls, n):
    """exact unnormalised integer moments (object array (24,)*n) and the normalisation prod_r |c_r|^2"""
    bound = 1
    for L in Ls: bound *= sum(abs(int(v)) for v in L) ** 2
    res = []; ps = []; prod = 1
    for p in PRIMES:
        res.append(moments_modp(Ls, n, p)); ps.append(p); prod *= p
        if prod > 4 * bound: break
    assert prod > 4 * bound, "not enough primes"
    flat = [r.ravel() for r in res]
    vals = np.empty(flat[0].size, dtype=object)
    for idx in range(flat[0].size):
        x = 0; M = 1
        for r, p in zip(flat, ps):
            t = ((int(r[idx]) - x) * pow(M, -1, p)) % p
            x = x + M * t; M *= p
        if x > M // 2: x -= M
        vals[idx] = x
    norm = 1
    for L in Ls: norm *= sum(int(v) * int(v) for v in L)
    return vals.reshape(res[0].shape), norm, len(ps)

def fourier_int(A):
    """integer F (23 x 24) and denominator dF: F/dF [(l,i,j), sigma] = f_l/(24 m_l) rho_l(sigma^-1)[i, j]"""
    rows = []
    for b in A.blocks:
        f = b['f']
        for i in range(f):
            for j in range(f):
                rows.append([fq(f, 24 * b['m']) * b['rho'][inv(s)][i, j] for s in S4])
    dF = 1
    for r in rows:
        for v in r: dF = dF * int(v.q) // math.gcd(dF, int(v.q))
    F = np.array([[int(v.p) * (dF // int(v.q)) for v in r] for r in rows], dtype=object)
    return F, dF

def hats_exact(m_int, norm, A, n):
    """site-symmetrised exact X_hat per sorted pattern (fmpq_mat), from integer moments"""
    F, dF = fourier_int(A)
    Z = m_int
    for s in range(n):
        Z = np.moveaxis(np.tensordot(F, Z, axes=([1], [s])), 0, s)
    Zs = None
    for perm in itertools.permutations(range(n)):
        Zp = Z.transpose(perm)
        Zs = Zp.copy() if Zs is None else Zs + Zp
    den = math.factorial(n) * dF ** n * norm
    nb = len(A.blocks); f = [b['f'] for b in A.blocks]
    offs = np.cumsum([0] + [x * x for x in f])
    hats = {}
    for lam in sorted_patterns(nb, n):
        sub = Zs[np.ix_(*[np.arange(offs[l], offs[l + 1]) for l in lam])]
        Y = pair_to_mat(sub[None], tuple(f[l] for l in lam))[0]
        D = Y.shape[0]
        hats[lam] = flint.fmpq_mat(D, D, [fq(int(v), den) for v in Y.ravel()])
    return hats

# ----------------------------------------------------------------------------- vectors
def cayley_int(rng):
    a, b, c = [int(z) for z in rng.integers(-3, 4, size=3)]
    S = flint.fmpq_mat(3, 3, [0, a, b, -a, 0, c, -b, -c, 0]); I = flint.fmpq_mat(3, 3, [1, 0, 0, 0, 1, 0, 0, 0, 1])
    Qm = (I - S) * (I + S).inv()
    d = 1
    for i in range(3):
        for j in range(3): d = d * int(Qm[i, j].q) // math.gcd(d, int(Qm[i, j].q))
    M = np.array([[int(Qm[i, j] * d) for j in range(3)] for i in range(3)], dtype=object)
    return M

def okr(*vs):
    out = np.array([1], dtype=object)
    for v in vs: out = np.kron(out, np.asarray(v, dtype=object))
    return out

PYTH = [(3, 4), (5, 12), (8, 15), (7, 24), (20, 21)]

def two_rows(rng):
    """an integer Cayley orthogonal matrix and an ordered pair of distinct row indices"""
    M = cayley_int(rng); i, j = [int(z) for z in rng.permutation(3)[:2]]
    return M[i], M[j]

def spinflip_Y(u):
    """v = Y u, Y|x> = (-1)^w(x) |xbar> on the qubit sector {0,1}^4 of (C^3)^(x)4 (u must vanish outside it)"""
    v = np.zeros(81, dtype=object)
    for x in itertools.product(range(3), repeat=4):
        idx = int(np.ravel_multi_index(x, (3,) * 4))
        if max(x) == 2:
            assert u[idx] == 0, "spin-flip vector leaves the qubit sector"; continue
        xb = tuple(1 - b for b in x); w = sum(x)
        v[int(np.ravel_multi_index(xb, (3,) * 4))] = (-1) ** w * u[idx]
    return v

def vectors(kind, rng):
    if kind == 'k2':
        a, a2 = two_rows(rng); b, b2 = two_rows(rng)
        while True:
            psi = rng.integers(-3, 4, size=9).astype(object)
            if flint.fmpz_mat(3, 3, [int(z) for z in psi]).det() != 0: break        # full Schmidt rank (entangled)
        return okr(a, b, psi), okr(a2, b2, psi)
    if kind == 'k4':
        rows = [two_rows(rng) for _ in range(4)]
        return okr(*[r[0] for r in rows]), okr(*[r[1] for r in rows])
    if kind == 'inner':
        a, a2 = two_rows(rng); rows = [two_rows(rng) for _ in range(3)]
        p = okr(*[r[0] for r in rows]); q = okr(*[r[1] for r in rows])
        r1, r2 = PYTH[int(rng.integers(0, len(PYTH)))]
        if rng.integers(0, 2): r1, r2 = r2, r1
        if rng.integers(0, 2): r2 = -r2
        phi = r1 * p + r2 * q; chi = -r2 * p + r1 * q
        return okr(a, phi), okr(a2, chi)
    if kind == 'sf':
        u = np.zeros(81, dtype=object); S = 0
        for x in itertools.product(range(2), repeat=4):
            if x[0] == 1 or x == (0, 0, 0, 0): continue                       # representatives of the 7 other pairs
            xb = tuple(1 - b for b in x)
            ux, uxb = [int(z) for z in rng.integers(-1, 2, size=2)]
            u[int(np.ravel_multi_index(x, (3,) * 4))] = ux; u[int(np.ravel_multi_index(xb, (3,) * 4))] = uxb
            S += (-1) ** sum(x) * ux * uxb
        u[0] = 1; u[int(np.ravel_multi_index((1, 1, 1, 1), (3,) * 4))] = -S
        v = spinflip_Y(u)
        assert sum(u * v) == 0          # u^T Y u = 2 sum over pairs (-1)^w(x) u_x u_xbar = 0
        return u, v
    if kind.startswith('near_'):
        # audit kind: a zero of kind <base>, scaled by 4, plus a random integer perturbation in {-1, 0, 1}^81 on u and
        # on v, then v orthogonalised against u exactly (small but nonzero q_4)
        u0, v0 = vectors(kind[5:], rng)
        p = 4 * u0 + rng.integers(-1, 2, size=81).astype(object); q = 4 * v0 + rng.integers(-1, 2, size=81).astype(object)
        return p, (p @ p) * q - (p @ q) * p
    if kind == 'ent':
        p = rng.integers(-2, 3, size=81).astype(object); q = rng.integers(-2, 3, size=81).astype(object)
        u = p; v = (p @ p) * q - (p @ q) * p
        return u, v
    raise ValueError(kind)

def qn_pm(u, v, t):
    """exact q_4(P_u + t P_v) from the 81 x 81 matrix"""
    Pm = np.outer(u, u) * mpq(1, int(sum(u * u))) + np.outer(v, v) * mpq(t, int(sum(v * v)))
    return XN.exact_qn(Pm, 4)

def main():
    kind, seed, name = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    t0 = time.time(); n = 4
    A, B1, A3 = site_algebras()
    vs = VarSpace(A, n, log=lambda s: None)
    rng = np.random.default_rng(seed)
    u, v = vectors(kind, rng)
    assert sum(u * v) == 0, "u, v not orthogonal"
    assert any(u) and any(v)
    if kind != 'ent' and not kind.startswith('near_'): assert sum(u * u) == sum(v * v)
    L = {'u': u, 'v': v}; hats = {}; npr = []
    for c in FAMS:
        m, norm, k = moments_int([L[ch] for ch in c], n); npr.append(k)
        hats[c] = hats_exact(m, norm, A, n)
        print(f"   {name}: family {c}: {k} primes; {time.time() - t0:.0f}s", flush=True)
    x, ok = vs.x_from_hats(hats)
    with open('lin4m.pkl', 'rb') as fh: lin = pickle.load(fh)
    xm = [mpq(int(q.p), int(q.q)) for q in x]
    Ex = [sum((row[i] * xm[i] for i in range(lin['N']) if row[i] != 0), mpq(0)) for row in lin['E']]
    eq_ok = all(Ex[i] == lin['e'][i] for i in range(len(Ex)))
    mom_ok = all(sum((r[i] * xm[i] for i in range(lin['N']) if r[i] != 0), mpq(0)) == rhs for r, rhs in zip(lin['mrows'], lin['mrhs']))
    obj = sum((lin['c'][i] * xm[i] for i in range(lin['N']) if lin['c'][i] != 0), mpq(0))
    objp = sum((lin['cplus'][i] * xm[i] for i in range(lin['N']) if lin['cplus'][i] != 0), mpq(0))
    # q_4 computed directly from the 81 x 81 matrices P = P_u - P_v (t = -1) and P_u + P_v (t = +1, check only)
    q4 = qn_pm(u, v, -1); q4p = qn_pm(u, v, 1)
    zero_ok = (obj == 0 and q4 == 0) if (kind != 'ent' and not kind.startswith('near_')) else None
    rec = dict(name=name, kind=kind, seed=seed, t=-1, u=[int(a) for a in u], v=[int(a) for a in v],
               x=[(int(q.p), int(q.q)) for q in x], in_space=ok, eq_ok=eq_ok, mom_ok=mom_ok,
               obj=(int(obj.numerator), int(obj.denominator)), q4=(int(q4.numerator), int(q4.denominator)),
               objplus=(int(objp.numerator), int(objp.denominator)), q4plus=(int(q4p.numerator), int(q4p.denominator)),
               zero_ok=zero_ok, primes=npr)
    os.makedirs('pts', exist_ok=True)
    with open(f'pts/{name}.pkl', 'wb') as fh: pickle.dump(rec, fh)
    xbits = max(max(int(abs(q.p)).bit_length(), int(q.q).bit_length()) for q in x)
    print(f"{name} ({kind}, seed {seed}): |u|^2 {int(sum(u * u))}, |v|^2 {int(sum(v * v))}, max|entry| "
          f"{max(max(abs(int(a)) for a in u), max(abs(int(a)) for a in v))}; in variable space {ok}; moments exact {mom_ok}; "
          f"equalities exact {eq_ok}; objective (t = -1) {float(obj):.12g} (exact {'0' if obj == 0 else 'nonzero'}); "
          f"q_4(P_u - P_v) direct {float(q4):.12g}; objective == q_4(P_u - P_v): {obj == q4}; t = +1 check: cplus.x == "
          f"q_4(P_u + P_v) = {float(q4p):.12g}: {objp == q4p}; x up to {xbits} bits; {time.time() - t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
