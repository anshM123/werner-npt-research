"""Exact removal of vacuous directions from the face blocks.

For a face block with data matrices F_i (symmetric, r x r), the common kernel K = cap_i ker F_i carries no constraint.
B (r' x r, integer, LLL-reduced rows) = basis of the row space of the stacked data [F_1; F_2; ...] = K^perp; the reduced
data F'_i = B F_i B^T define the constraint B F(x) B^T >= 0, implied by F(x) >= 0 (valid for any B), and a dual S' >= 0
of the reduced block is the dual S = B^T S' B >= 0 of the original one.  Strict dual feasibility is preserved
(S' = M S M^T with M = (B B^T)^{-1} B).
Rows of the stacked matrix independent modulo p = 1048573 are selected; completeness of the row space is then checked
EXACTLY: F_i Kb = 0 for every data matrix, Kb an exact basis of ker B.
usage: python span4.py   (reads facedata4.pkl, writes facedata4r.pkl)"""
import sys, time, pickle, math, numpy as np, flint
from gmpy2 import mpq
from build_lin4 import independent_rows_modp, P_SEL

def stacked_int(b):
    """list of integer matrices (object arrays r x r), one per data matrix (rows scaled by positive integers)"""
    r = b['r']; out = []
    if b['kind'] == 'COV':
        for gi, idx, Fh, den in b['data']:
            for k in range(Fh.shape[0]): out.append(Fh[k])
    else:
        gi, idx, mats = b['data'][0]
        for M in mats:
            d = 1
            for row in M:
                for (a, q) in row: d = d * q // math.gcd(d, q)
            out.append(np.array([[a * (d // q) for (a, q) in row] for row in M], dtype=object))
    return out

def main():
    t0 = time.time()
    fd = pickle.load(open('facedata4.pkl', 'rb'))
    newblocks = []; tot = [0, 0]; nred = 0
    for b in fd['blocks']:
        r = b['r']; mats = stacked_int(b)
        St = np.vstack(mats)                                   # (sum K r) x r object ints
        Stp = np.array([[int(v) % P_SEL for v in row] for row in St], dtype=np.int64)
        sel = independent_rows_modp(Stp)
        rows = St[sel]
        Rq = flint.fmpq_mat(len(sel), r, [flint.fmpq(int(v)) for v in rows.ravel()])
        RR, rk = Rq.rref()
        assert rk == len(sel)
        d = 1
        for i in range(rk):
            for j in range(r): d = d * int(RR[i, j].q) // math.gcd(d, int(RR[i, j].q))
        Z = flint.fmpz_mat(rk, r, [int((RR[i, j] * d).p) for i in range(rk) for j in range(r)]).lll()
        B = np.array([[int(Z[i, j]) for j in range(r)] for i in range(rk)], dtype=object)
        # exact completeness check: F_i Kb = 0 for an exact basis Kb of ker B
        Kb, nul = flint.fmpz_mat(rk, r, [int(v) for v in B.ravel()]).nullspace()
        Kbo = np.array([[int(Kb[i, j]) for j in range(nul)] for i in range(r)], dtype=object).reshape(r, nul)
        if nul:
            prod = St.dot(Kbo)
            assert all(v == 0 for v in prod.ravel()), f"{b['name']}: row space incomplete modulo p"
        tot[0] += r; tot[1] += rk
        if rk < r: nred += 1
        if rk == 0: continue
        nb = dict(name=b['name'], kind=b['kind'], r=rk, B=B)
        if b['kind'] == 'COV':
            nb['data'] = [(gi, idx, np.array([B.dot(Fh[k]).dot(B.T) for k in range(Fh.shape[0])], dtype=object), den)
                          for gi, idx, Fh, den in b['data']]
        else:
            gi, idx, mats0 = b['data'][0]
            newm = []
            for M in mats0:
                Mq = flint.fmpq_mat(r, r, [flint.fmpq(a, q) for row in M for (a, q) in row])
                Bq = flint.fmpq_mat(rk, r, [int(v) for v in B.ravel()])
                F = Bq * Mq * Bq.transpose()
                newm.append([[(int(F[i, j].p), int(F[i, j].q)) for j in range(rk)] for i in range(rk)])
            nb['data'] = [(gi, idx, newm)]
        newblocks.append(nb)
    bbits = max(max(abs(int(v)) for v in b['B'].ravel()).bit_length() for b in newblocks)
    nbits = max(max(abs(int(v)).bit_length() for gi, idx, Fh, den in b['data'] for v in Fh.ravel()) for b in newblocks if b['kind'] == 'COV')
    pickle.dump(dict(blocks=newblocks, primes=fd['primes']), open('facedata4r.pkl', 'wb'))
    print(f"[span] {len(fd['blocks'])} face blocks, {nred} with vacuous directions; total size {tot[0]} -> {tot[1]}; "
          f"{len(newblocks)} blocks kept (largest {max(b['r'] for b in newblocks)}, svec {sum(b['r'] * (b['r'] + 1) // 2 for b in newblocks)}); "
          f"B entries <= {bbits} bits, COV numerators <= {nbits} bits; row-space completeness checked exactly; {time.time() - t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
