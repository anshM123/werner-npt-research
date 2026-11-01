"""Which lifted block directions are killed by the equalities (diagnostic; not part of the certificate).

For every lifted block k and every diagonal entry a of its (+,+) part (functional f(x) = F+'_k(x)[a, a]) and of its
(-,-) part (f(x) = F-'_k(x)[a, a]), test whether f vanishes on the affine space E x = e, i.e. f = lambda^T E with
lambda . e = 0.  E (linL.pkl) is in reduced row echelon form, so lambda_j = f[pivot_j] and the test is
f - sum_j lambda_j E_j = 0 and lambda . e = 0.  Screening modulo the prime 134217689 (a nonzero residue proves that f
does not vanish there); every candidate is then confirmed in exact rational arithmetic.  For such a direction the
whole row of the lifted block vanishes at every feasible point (PSD), so the corresponding dual entries are not pinned
down by the identities at true points (they are by the exact identities for all x).
usage: python affineL.py"""
import pickle, time, numpy as np
from gmpy2 import mpq

P = 134217689

def main():
    t0 = time.time()
    lin = pickle.load(open('linL.pkl', 'rb')); fd = pickle.load(open('facedataLr.pkl', 'rb'))
    N = lin['N']; E = lin['E']; e = lin['e']; neq = len(E); blocks = fd['blocks']
    piv = []
    for j in range(neq):
        p = next(i for i in range(N) if E[j][i] != 0); assert E[j][p] == 1 and all(E[jj][p] == 0 for jj in range(neq) if jj != j)
        piv.append(p)
    md = lambda q: (int(q.numerator) % P) * pow(int(q.denominator) % P, -1, P) % P
    Ep = np.array([[md(v) for v in row] for row in E], dtype=np.int64); ep = np.array([md(v) for v in e], dtype=np.int64)
    cand = []; total = 0; per_block = {}
    for k, b in enumerate(blocks):
        for part, r in (('+', b['rp']), ('-', b['rm'])):
            for a in range(r):
                total += 1
                f = np.zeros(N, dtype=np.int64)
                for gi, idx, Fp, Fm, G, D in b['data']:
                    M = Fp if part == '+' else Fm
                    if M.size == 0: continue
                    inv = pow(int(D) % P, -1, P)
                    for kk, i in enumerate(idx):
                        v = int(M[kk, a, a])
                        if v: f[i] = (f[i] + (v % P) * inv) % P
                lam = f[piv]
                res = (f - (lam[:, None] * Ep % P).sum(axis=0)) % P
                if not res.any() and int((lam * ep % P).sum() % P) == 0: cand.append((k, part, a))
    confirmed = []
    for (k, part, a) in cand:
        b = blocks[k]; f = [mpq(0)] * N
        for gi, idx, Fp, Fm, G, D in b['data']:
            M = Fp if part == '+' else Fm
            if M.size == 0: continue
            for kk, i in enumerate(idx):
                if M[kk, a, a]: f[i] += mpq(int(M[kk, a, a]), D)
        lam = [f[p] for p in piv]
        ok = all(f[i] - sum((lam[j] * E[j][i] for j in range(neq) if lam[j] != 0), mpq(0)) == 0 for i in range(N)) and \
             sum((lam[j] * e[j] for j in range(neq)), mpq(0)) == 0
        if ok: confirmed.append((k, part, a)); per_block.setdefault(b['name'], []).append(f"{part}{a}")
    print(f"[affine] {total} diagonal functionals of the lifted blocks; {len(cand)} vanish on E x = e modulo {P}, "
          f"{len(confirmed)} confirmed exactly, in {len(per_block)} blocks; {time.time() - t0:.0f}s", flush=True)
    for name, lst in sorted(per_block.items()): print(f"   {name}: {' '.join(lst)}")
    pickle.dump(dict(confirmed=confirmed), open('affineL.pkl', 'wb'))

if __name__ == "__main__":
    main()
