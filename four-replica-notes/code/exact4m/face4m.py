"""Face of the n = 4, t = -1 relaxation at the exact zeros of q_4(P_u - P_v) (stage 1 of the covariance data).

x* = sum_z lambda_z x_z with lambda_z = denominator of x_z (positive weights), an integer vector a.  Since every
W(x_z) is PSD (true points), ker W(x*) = intersection of ker W(x_z).
  X blocks:   Y_g(x*) = sum_k a_k Y_k exactly (flint), kernel exact.
  COV blocks: for every compressed block (xy, mu, sigma), the numerator
                 What(x*) = sym( [[ U^T Gt Yhat_rs(x*) U ]]_{r,s} ),   Yhat_rs(x*) = sum_g (LX/dX_g) sum_k a_k Ytil_k
              (Ytil_k = exact integer tree image of parameter k, float64-exact) is computed modulo 27-bit primes and
              recovered by CRT with an a-priori bound; W(x*) = What / (2 dT^n LX) (positive scale).
Kernels: exact integer bases (flint nullspace), LLL-reduced.  Saves face4m.pkl (or the file given by -o; with
--nosave only the dimensions are reported, e.g. for saturation checks with subsets of the zeros).
usage: python face4m.py [-o out.pkl | --nosave] zero_name [zero_name ...]"""
import sys, os, time, pickle, math, itertools, numpy as np, flint
from gmpy2 import mpq
from core4 import *
from engine4 import *
from symm4 import cov_classes
from points4m import PRIMES

def load_x(name):
    with open(f'pts/{name}.pkl', 'rb') as fh: rec = pickle.load(fh)
    assert rec['in_space'] and rec['eq_ok'] and rec['mom_ok'] and rec['obj'] == (0, 1), f"{name} is not a checked zero"
    assert rec.get('t') == -1 and rec['kind'] != 'ent' and rec['q4'] == (0, 1) and rec['zero_ok'], \
        f"{name}: not a checked zero of q_4(P_u - P_v)"
    return rec

def int_kernel(Wq):
    """integer LLL-reduced basis (D x r object array) of the kernel of an fmpz_mat / fmpq_mat (square)"""
    if isinstance(Wq, flint.fmpq_mat):
        d = 1
        for i in range(Wq.nrows()):
            for j in range(Wq.ncols()): d = d * int(Wq[i, j].q) // math.gcd(d, int(Wq[i, j].q))
        Wq = flint.fmpz_mat(Wq.nrows(), Wq.ncols(), [int((Wq[i, j] * d).p) for i in range(Wq.nrows()) for j in range(Wq.ncols())])
    X, nul = Wq.nullspace()
    D = Wq.ncols()
    if nul == 0: return np.zeros((D, 0), dtype=object)
    rows = []
    for k in range(nul):
        v = [int(X[i, k]) for i in range(D)]
        g = 0
        for a in v: g = math.gcd(g, a)
        rows.append([a // g for a in v])
    Z = flint.fmpz_mat(nul, D, [a for r in rows for a in r]).lll()
    N = np.array([[int(Z[k, i]) for k in range(nul)] for i in range(D)], dtype=object)
    # exactness: W N = 0
    Nq = flint.fmpz_mat(D, nul, [int(v) for v in N.ravel()])
    assert (Wq * Nq).is_zero(), 'kernel basis check failed'
    return N

def main():
    args = sys.argv[1:]; out = 'face4m.pkl'
    if args and args[0] == '-o': out = args[1]; args = args[2:]
    elif args and args[0] == '--nosave': out = None; args = args[1:]
    names = args
    t0 = time.time(); n = 4
    A, B1, A3 = site_algebras()
    vs = VarSpace(A, n, log=lambda s: None); N = vs.N
    recs = [load_x(nm) for nm in names]
    a = [0] * N
    for rec in recs:
        lam_z = 1
        for (p_, q_) in rec['x']: lam_z = lam_z * q_ // math.gcd(lam_z, q_)
        for i, (p_, q_) in enumerate(rec['x']): a[i] += p_ * (lam_z // q_)
    abits = max(abs(v).bit_length() for v in a)
    print(f"[face] {len(recs)} zeros: {', '.join(names)}; x* integer, up to {abits} bits; {time.time() - t0:.1f}s", flush=True)
    # ---------------- X blocks
    xker = {}; dims_x = [0, 0]
    for key in vs.order:
        g = vs.groups[key]; t = g['t']
        Y = flint.fmpq_mat(t, t)
        for k, i in enumerate(g['idx']):
            if a[i] != 0: Y = Y + g['Ys'][k] * a[i]
        Nk = int_kernel(Y)
        xker[key] = Nk; dims_x[0] += t; dims_x[1] += Nk.shape[1]
    print(f"[face] X blocks: dimension {dims_x[0]} -> {dims_x[1]}; {time.time() - t0:.1f}s", flush=True)
    # ---------------- COV blocks
    LX = 1
    for key in vs.order: LX = LX * vs.groups[key]['dX'] // math.gcd(LX, vs.groups[key]['dX'])
    entries = cov_entries()
    maps = {}
    for (xy, r, s_, c, P) in entries:
        Pm = PERM[P]
        Ti, dT = site_map_int(A, B1, lambda X, Pm=Pm: realign1(Pm @ X @ Pm.T))
        maps[(xy, r, s_)] = (c, Ti, dT)
    dTs = set(m[2] for m in maps.values()); assert len(dTs) == 1, "covariance site maps with different denominators"
    dT = dTs.pop(); tau = max(rowsum_bound(m[1]) for m in maps.values())
    mus = sorted_patterns(len(B1.blocks), n)
    classes = {(xy, mu): cov_classes(B1, xy, mu) for xy in ['uu', 'uv', 'vv'] for mu in mus}
    Gts = {mu: np.array(fm_to_int(fm_kron([B1.blocks[q]['G'] for q in mu]), 1), dtype=object) for mu in mus}
    # a-priori bound for |Yhat(x*)| entries and for the compressed numerators
    ybound = 0
    for key in vs.order:
        g = vs.groups[key]
        amax = max(abs(a[i]) for i in g['idx'])
        ybound += (LX // g['dX']) * amax * g['K'] * len(g['orbit']) * tau ** n * g['Xmax']
    cbound = 0
    for (xy, mu), cl in classes.items():
        for lab, dsig, U in cl:
            GU = Gts[mu].dot(U)
            cbound = max(cbound, int(np.abs(GU).sum(axis=0).max()) * int(np.abs(U).sum(axis=0).max()))
    wbound = 2 * 2 * cbound * ybound          # sym sum of two terms, (r,s) entries each bounded by ybound
    primes = []; prod = 1
    for p in PRIMES:
        primes.append(p); prod *= p
        if prod > 2 * wbound + 1: break
    assert prod > 2 * wbound + 1, "not enough primes for W(x*)"
    print(f"[face] COV: {sum(len(v) for v in classes.values())} class blocks; |Yhat(x*)| <= 2^{ybound.bit_length()}, "
          f"|What| <= 2^{wbound.bit_length()}: {len(primes)} primes; {time.time() - t0:.1f}s", flush=True)
    # accumulators: Yhat_rs(x*) mod p per (xy, mu, r, s)
    acc = {}
    for (xy, r, s_), (c, Ti, _) in maps.items():
        for key in vs.order:
            if key[0] != c: continue
            g = vs.groups[key]
            level = tree(vs.pair_float(key), key[1], Ti, len(B1.blocks))
            bnd = g['Xmax'] * tau ** n * len(g['orbit'])
            ak = [a[i] * (LX // g['dX']) for i in g['idx']]
            if not any(ak): continue
            for mu in mus:
                Y = assemble(level, g['orbit'], mu)
                if Y is None: continue
                check_exact(Y, bnd)
                gs = tuple(B1.blocks[q]['f'] for q in mu)
                Yi = pair_to_mat(Y, gs).astype(np.int64)                  # (K, Dm, Dm) exact
                for p in primes:
                    akp = np.array([v % p for v in ak], dtype=np.int64)
                    S = np.tensordot(akp, Yi % p, axes=1) % p          # K <= 96 terms < 2^54 each
                    kk = (xy, mu, r, s_, p)
                    acc[kk] = (acc[kk] + S) % p if kk in acc else S
        print(f"   entry {xy}{'uv'[r]}{'uv'[s_]} (family {c}); {time.time() - t0:.1f}s", flush=True)
    # compress, CRT, kernels
    cker = {}; dims_c = [0, 0]; nzero = 0
    for (xy, mu), cl in classes.items():
        Dm = Gts[mu].shape[0]
        for lab, dsig, U in cl:
            m = U.shape[1]
            res = []
            for p in primes:
                Up = np.array([[int(v) % p for v in row] for row in U], dtype=np.int64)
                GUp = np.array([[int(v) % p for v in row] for row in Gts[mu].dot(U)], dtype=np.int64)
                Wt = np.zeros((2 * m, 2 * m), dtype=np.int64)
                for r in range(2):
                    for s_ in range(2):
                        S = acc.get((xy, mu, r, s_, p))
                        if S is None: continue
                        Wt[r * m:(r + 1) * m, s_ * m:(s_ + 1) * m] = ((GUp.T @ S) % p @ Up) % p
                res.append((Wt + Wt.T) % p)
            # CRT to symmetric residues
            Wex = np.empty((2 * m, 2 * m), dtype=object)
            for i in range(2 * m):
                for j in range(2 * m):
                    x = 0; M = 1
                    for rr, p in zip(res, primes):
                        t = ((int(rr[i, j]) - x) * pow(M, -1, p)) % p
                        x += M * t; M *= p
                    if x > M // 2: x -= M
                    Wex[i, j] = x
            if all(v == 0 for v in Wex.ravel()): nzero += 1
            Wz = flint.fmpz_mat(2 * m, 2 * m, [int(v) for v in Wex.ravel()])
            Nk = int_kernel(Wz)
            # PSD sanity (float) of W(x*)
            ev = np.linalg.eigvalsh(np.array(Wex, dtype=float))
            assert ev.min() > -1e-9 * max(1.0, abs(ev).max()), f"W(x*) not PSD for {xy}{mu}{lab}: {ev.min()}"
            cker[(xy, mu, lab)] = Nk; dims_c[0] += 2 * m; dims_c[1] += Nk.shape[1]
    nbits = max((max(abs(int(v)) for v in Nk.ravel()).bit_length() for Nk in cker.values() if Nk.size), default=0)
    print(f"[face] COV blocks: dimension {dims_c[0]} -> {dims_c[1]} ({nzero} blocks with W(x*) = 0); kernel bases up to "
          f"{nbits} bits; {time.time() - t0:.1f}s", flush=True)
    if out is None:
        print(f"[face] not saved (--nosave); {time.time() - t0:.1f}s", flush=True); return
    with open(out, 'wb') as fh:
        pickle.dump(dict(names=names, a=a, xker=xker, cker=cker, LX=LX, dT=dT,
                         U={k: [lab_U for lab_U in v] for k, v in classes.items()}), fh)
    print(f"saved {out}; {time.time() - t0:.1f}s", flush=True)

if __name__ == "__main__":
    main()
