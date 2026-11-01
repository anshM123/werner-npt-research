"""Exact points of the relaxation at n = 4: zeros of q_4 (t = +1) and an entangled test point.

  vectors   integer u, v in Z^81 with <u, v> = 0 and |u| = |v| (zeros) or u = p, v = (p.p) q - (p.q) p (test point);
  moments   m[sigma] = Tr[(x)_r |c_r><c_r| (sigma_1 (x) ... (x) sigma_4)] (integer, unnormalised) modulo 27-bit
            primes (pairwise einsum with an asserted bound: summed terms * p^2 < 2^63) and CRT with the bound
            prod_r ||c_r||_1^2 (as in ../relaxations/exact_coreN.moments_exact);
  hats      exact integer Fourier transform per site (F[(l,i,j), sigma] = f_l/(24 m_l) rho_l(sigma^-1)[i,j]),
            site symmetrisation, exact X_hat per sorted pattern;
  x         exact projection onto the parameters with exactness checks (VarSpace.x_from_hats).
Saves pts/<name>.pkl with x (list of (num, den)), the vectors, and the exact checks.
usage: python points4.py <kind> <seed> <name>     kind in k1, k3, rk1, rk3, rint, ent"""
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

def vectors(kind, rng):
    def nz(v):
        if not any(v): v[0] = 1
        return v
    if kind in ('k1', 'rk1'):
        M = cayley_int(rng); psi = nz(rng.integers(-3, 4, size=27).astype(object))
        u, v = okr(M[0], psi), okr(M[1], psi)
    elif kind in ('k3', 'rk3'):
        M1, M2, M3 = cayley_int(rng), cayley_int(rng), cayley_int(rng); psi = nz(rng.integers(-3, 4, size=3).astype(object))
        u, v = okr(M1[0], M2[0], M3[0], psi), okr(M1[1], M2[1], M3[1], psi)
    elif kind == 'rint':
        Ma, Mb, Mc = cayley_int(rng), cayley_int(rng), cayley_int(rng); psi = nz(rng.integers(-3, 4, size=3).astype(object))
        wp = okr(Mb[0], Mc[0], psi); wm = okr(Mb[1], Mc[1], psi)
        p_, q_ = PYTH[int(rng.integers(0, len(PYTH)))]
        if rng.integers(0, 2): p_, q_ = q_, p_
        s = 1 if rng.integers(0, 2) else -1
        phi = p_ * wp + q_ * s * wm; chi = q_ * wp + p_ * s * wm
        u, v = okr(Ma[0], phi), okr(Ma[1], chi)
    elif kind == 'ent':
        p = rng.integers(-2, 3, size=81).astype(object); q = rng.integers(-2, 3, size=81).astype(object)
        u = p; v = (p @ p) * q - (p @ q) * p
        return u, v
    else:
        raise ValueError(kind)
    if kind.startswith('r'):
        r1, r2 = PYTH[int(rng.integers(0, len(PYTH)))]
        u, v = r1 * u + r2 * v, -r2 * u + r1 * v
    return u, v

def main():
    kind, seed, name = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    t0 = time.time(); n = 4
    A, B1, A3 = site_algebras()
    vs = VarSpace(A, n, log=lambda s: None)
    rng = np.random.default_rng(seed)
    u, v = vectors(kind, rng)
    assert sum(u * v) == 0, "u, v not orthogonal"
    if kind != 'ent': assert sum(u * u) == sum(v * v)
    L = {'u': u, 'v': v}; hats = {}; npr = []
    for c in FAMS:
        m, norm, k = moments_int([L[ch] for ch in c], n); npr.append(k)
        hats[c] = hats_exact(m, norm, A, n)
        print(f"   {name}: family {c}: {k} primes; {time.time() - t0:.0f}s", flush=True)
    x, ok = vs.x_from_hats(hats)
    with open('lin4.pkl', 'rb') as fh: lin = pickle.load(fh)
    xm = [mpq(int(q.p), int(q.q)) for q in x]
    Ex = [sum((row[i] * xm[i] for i in range(lin['N']) if row[i] != 0), mpq(0)) for row in lin['E']]
    eq_ok = all(Ex[i] == lin['e'][i] for i in range(len(Ex)))
    mom_ok = all(sum((r[i] * xm[i] for i in range(lin['N']) if r[i] != 0), mpq(0)) == rhs for r, rhs in zip(lin['mrows'], lin['mrhs']))
    obj = sum((lin['c'][i] * xm[i] for i in range(lin['N']) if lin['c'][i] != 0), mpq(0))
    # q_4(P) computed directly from P
    Pm = np.outer(u, u) * mpq(1, int(sum(u * u))) + np.outer(v, v) * mpq(1, int(sum(v * v)))
    q4 = XN.exact_qn(Pm, n)
    rec = dict(name=name, kind=kind, seed=seed, u=[int(a) for a in u], v=[int(a) for a in v],
               x=[(int(q.p), int(q.q)) for q in x], in_space=ok, eq_ok=eq_ok, mom_ok=mom_ok,
               obj=(int(obj.numerator), int(obj.denominator)), q4=(int(q4.numerator), int(q4.denominator)), primes=npr)
    os.makedirs('pts', exist_ok=True)
    with open(f'pts/{name}.pkl', 'wb') as fh: pickle.dump(rec, fh)
    xbits = max(max(int(abs(q.p)).bit_length(), int(q.q).bit_length()) for q in x)
    print(f"{name} ({kind}, seed {seed}): in variable space {ok}; moments exact {mom_ok}; equalities exact {eq_ok}; "
          f"objective {float(obj):.12g} (exact {'0' if obj == 0 else 'nonzero'}); q_4(P) direct {float(q4):.12g}; "
          f"objective == q_4(P): {obj == q4}; x up to {xbits} bits; {time.time() - t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
