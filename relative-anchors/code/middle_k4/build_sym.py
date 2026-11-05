"""Exact symbolic construction of the letter-1 and type I torus blocks of (R_inf)_1 and (R_inf)_2 at k = 4.

ylib.Builder over Q(s1, s2, s3) with s0 = 1 (rescaled basis; entries homogeneous of degree n1 = j + 1 - |T| in
(s0, ..., s3), so s0 = 1 loses nothing).  For each block the script checks exactly that the known kernel vectors are
annihilated (as identities of rational functions), and stores the full block and the gauge-fixed indices in
blocks_sym.pkl.  Also prints the common denominator of each block and its factorisation.
Usage: python build_sym.py"""
import time, pickle
import flint
from k4lib import block_states, kernel_vectors, gauge_keep, sym_builder, BLOCKS, label

t0 = time.time()
B, K = sym_builder(dehom=True)
s = [K(1)] + list(K.gens)
out = {}
ok = True
for (j, T) in BLOCKS:
    S = block_states(j, T)
    M = B.matrix(S)
    n = len(S)
    sym = all(M[a][b] == M[b][a] for a in range(n) for b in range(n))
    kv = kernel_vectors(S, j, T, s)
    kz = all(all(sum((M[i][l] * v[l] for l in range(n) if v[l] != 0), K(0)) == 0 for i in range(n)) for v in kv)
    keep = gauge_keep(S, j, T)
    den = None
    for a in range(n):
        for b in range(n):
            d = M[a][b].d
            if den is None:
                den = d
            else:
                g = den.gcd(d)
                den = den * (d / g)
    c, fac = den.factor()
    facs = ' * '.join(f'[deg {f.total_degree()}, {len(f)} t, pos {all(x > 0 for x in f.coeffs())}]^{e}' for f, e in fac)
    maxdeg = max(M[a][b].n.total_degree() for a in range(n) for b in range(n))
    print(f'j = {j}, T = {T}: size {n}, symmetric {sym}, known kernel {len(kv)} annihilated exactly: {kz}; '
          f'gauge-fixed size {len(keep)}; common den = {c} * {facs}; max num degree {maxdeg}  [{time.time() - t0:.1f}s]',
          flush=True)
    ok &= sym and kz and len(keep) == n - len(kv)
    out[(j, T)] = dict(states=S, labels=[label(x) for x in S], keep=keep,
                       M=[[(M[a][b].n.to_dict(), M[a][b].d.to_dict()) for b in range(n)] for a in range(n)])
pickle.dump(dict(names=['s1', 's2', 's3'], blocks=out), open('blocks_sym.pkl', 'wb'))
print(f'time {time.time() - t0:.1f}s')
print('PASS' if ok else 'FAIL')
