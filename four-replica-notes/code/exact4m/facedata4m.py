"""Exact face-reduced constraint data at n = 4, t = -1 (stage 2; reads face4m.pkl).

COV block k = (xy, mu, sigma) with kernel basis N = [N_u; N_v] (2m x r, integer):  for parameter i of group g,
    F_k(e_i) = Fhat_{k,i} / (2 dT^n dX_g),   Fhat = F0 + F0^T,   F0 = sum_{r,s} (Gt U N_r)^T Ytil_rs,i (U N_s),
Ytil_rs,i the exact integer tree image of parameter i in the sorted target mu (float64-exact).  F0 is computed modulo
27-bit primes; Fhat is recovered by CRT with the a-priori bound
    |Fhat| <= 2 sum_{r,s} max|Ytil_rs| colsum(Gt U N_r) colsum(U N_s),   max|Ytil| <= |orbit| tau^n Xmax_g.
X block of group g with kernel N_g (t x r_g):  F_g(e_k) = N_g^T Y_k N_g (exact, flint).
Checks: F_k(x_z) = 0 exactly for every zero z (the kernel bases lie in the face), and F(x*) = 0.
Saves facedata4m.pkl:  blocks = list of dict(name, r, kind, data) with data = list of (group index, idx list,
numerator object array (K, r, r), denominator) for COV, and (idx, list of fmpq-as-(num,den) r x r) for X blocks.
usage: python facedata4m.py"""
import sys, os, time, pickle, math, itertools, numpy as np, flint
from gmpy2 import mpq
from core4 import *
from engine4 import *
from points4m import PRIMES

def crt_array(res, primes):
    """CRT of a list of int64 residue arrays (same shape) to symmetric representatives (object array)"""
    shape = res[0].shape; flat = [r.ravel() for r in res]
    M = 1; x = np.zeros(flat[0].size, dtype=object)
    for r, p in zip(flat, primes):
        inv_ = pow(M % p, -1, p)
        t = ((r.astype(object) - x) % p) * inv_ % p
        x = x + M * t; M *= p
    half = M // 2
    x = np.where(x > half, x - M, x)
    return x.reshape(shape)

def main():
    t0 = time.time(); n = 4
    A, B1, A3 = site_algebras()
    vs = VarSpace(A, n, log=lambda s: None); N = vs.N
    with open('face4m.pkl', 'rb') as fh: face = pickle.load(fh)
    xker, cker, Ub = face['xker'], face['cker'], face['U']
    mus = sorted_patterns(len(B1.blocks), n)
    Gts = {mu: np.array(fm_to_int(fm_kron([B1.blocks[q]['G'] for q in mu]), 1), dtype=object) for mu in mus}
    maps = {}
    for (xy, r, s_, c, P) in cov_entries():
        Pm = PERM[P]
        Ti, dT = site_map_int(A, B1, lambda X, Pm=Pm: realign1(Pm @ X @ Pm.T))
        maps[(xy, r, s_)] = (c, Ti, dT)
    dT = face['dT']; assert all(m[2] == dT for m in maps.values())
    tau = max(rowsum_bound(m[1]) for m in maps.values())
    # kept COV blocks and their compression matrices
    kept = []          # (xy, mu, lab, m, r, L[0..1], R[0..1])
    for (xy, mu), cl in Ub.items():
        for lab, dsig, U in cl:
            Nk = cker[(xy, mu, lab)]; m = U.shape[1]
            if Nk.shape[1] == 0: continue
            NL = [Nk[:m], Nk[m:]]
            L = [Gts[mu].dot(U).dot(NL[r]) for r in range(2)]          # Gt U N_r: Dm x r
            R = [U.dot(NL[s_]) for s_ in range(2)]
            kept.append(dict(xy=xy, mu=mu, lab=lab, m=m, r=Nk.shape[1], L=L, R=R))
    kidx = {(b['xy'], b['mu']): [] for b in kept}
    for bi, b in enumerate(kept): kidx[(b['xy'], b['mu'])].append(bi)
    # CRT bound over all blocks and groups
    gmax = max(len(vs.groups[k]['orbit']) * vs.groups[k]['Xmax'] for k in vs.order) * tau ** n
    fb = 0
    for b in kept:
        cs = lambda M: int(np.abs(M).sum(axis=0).max()) if M.size else 0
        tot = sum(cs(b['L'][r]) * cs(b['R'][s_]) for r in range(2) for s_ in range(2))
        fb = max(fb, 2 * gmax * tot)
    primes = []; prod = 1
    for p in PRIMES:
        primes.append(p); prod *= p
        if prod > 2 * fb + 1: break
    assert prod > 2 * fb + 1
    nsvec = sum(b['r'] * (b['r'] + 1) // 2 for b in kept)
    print(f"[data] {len(kept)} COV face blocks (max r {max(b['r'] for b in kept)}, svec {nsvec}); |Fhat| <= 2^{fb.bit_length()}: "
          f"{len(primes)} primes; {time.time() - t0:.1f}s", flush=True)
    # stacked compression matrices per (xy, mu): Lall_r = [Gt U N_r of every kept block] (Dm x sum r)
    def modp(M, p): return np.array([[int(v) % p for v in row] for row in M], dtype=np.int64).reshape(M.shape)
    Lst = {}; Rst = {}; offs = {}
    for (xy, mu), bis in kidx.items():
        o = [0]
        for bi in bis: o.append(o[-1] + kept[bi]['r'])
        offs[(xy, mu)] = o
        for p in primes:
            for r in range(2):
                Lst[(xy, mu, r, p)] = modp(np.hstack([kept[bi]['L'][r] for bi in bis]), p)
                for bi in bis: Rst[(bi, r, p)] = modp(kept[bi]['R'][r], p)
    acc = {}       # (bi, group key) -> dict prime -> int64 (K, r, r)
    for (xy, r, s_), (c, Ti, _) in maps.items():
        for key in vs.order:
            if key[0] != c: continue
            g = vs.groups[key]
            level = tree(vs.pair_float(key), key[1], Ti, len(B1.blocks))
            bnd = g['Xmax'] * tau ** n * len(g['orbit'])
            for mu in mus:
                if (xy, mu) not in kidx: continue
                Y = assemble(level, g['orbit'], mu)
                if Y is None: continue
                check_exact(Y, bnd)
                gs = tuple(B1.blocks[q]['f'] for q in mu)
                Yi = pair_to_mat(Y, gs).astype(np.int64)
                o = offs[(xy, mu)]
                for p in primes:
                    Yp = Yi % p
                    T1 = np.matmul(Lst[(xy, mu, r, p)].T, Yp) % p            # (K, sum r, Dm); sums of <= 81 terms < 2^54
                    for j, bi in enumerate(kidx[(xy, mu)]):
                        T2 = np.matmul(T1[:, o[j]:o[j + 1], :], Rst[(bi, s_, p)]) % p
                        kk = (bi, key)
                        if kk not in acc: acc[kk] = {q: np.zeros((g['K'], kept[bi]['r'], kept[bi]['r']), dtype=np.int64) for q in primes}
                        acc[kk][p] = (acc[kk][p] + T2) % p
        print(f"   entry {xy}{'uv'[r]}{'uv'[s_]} (family {c}); {time.time() - t0:.1f}s", flush=True)
    # symmetrise (Fhat = F0 + F0^T) and CRT
    blocks = []
    for bi, b in enumerate(kept):
        data = []
        for key in vs.order:
            kk = (bi, key)
            if kk not in acc: continue
            res = [(acc[kk][p] + acc[kk][p].transpose(0, 2, 1)) % p for p in primes]
            Fh = crt_array(res, primes)
            if all(v == 0 for v in Fh.ravel()): continue
            g = vs.groups[key]
            data.append((vs.order.index(key), list(g['idx']), Fh, 2 * dT ** n * g['dX']))
        if data:
            blocks.append(dict(name=f"COV{b['xy']}{b['mu']}{b['lab']}", kind='COV', r=b['r'], data=data))
    del acc
    print(f"[data] COV: {len(blocks)} nonzero face blocks after CRT; {time.time() - t0:.1f}s", flush=True)
    # X blocks
    nx = 0
    for key in vs.order:
        Ng = xker[key]
        if Ng.shape[1] == 0: continue
        g = vs.groups[key]; t = g['t']; rg = Ng.shape[1]
        Nq = flint.fmpq_mat(t, rg, [int(v) for v in Ng.ravel()])
        mats = []
        for k in range(g['K']):
            F = Nq.transpose() * g['Ys'][k] * Nq
            mats.append([[(int(F[a, b_].p), int(F[a, b_].q)) for b_ in range(rg)] for a in range(rg)])
        if all(v[0] == 0 for M in mats for row in M for v in row): continue
        blocks.append(dict(name=f"X{key[0]}{key[1]}", kind='X', r=rg, data=[(vs.order.index(key), list(g['idx']), mats)]))
        nx += 1
    print(f"[data] X: {nx} nonzero face blocks; {time.time() - t0:.1f}s", flush=True)
    # exact face checks at the zeros
    zx = []
    for nm in face['names']:
        with open(f'pts/{nm}.pkl', 'rb') as fh: rec = pickle.load(fh)
        zx.append((nm, [mpq(p_, q_) for (p_, q_) in rec['x']]))
    bad = 0
    for nm, x in zx:
        for blk in blocks:
            r = blk['r']; tot = np.full((r, r), mpq(0), dtype=object)
            if blk['kind'] == 'COV':
                for gi, idx, Fh, den in blk['data']:
                    xv = np.array([x[i] for i in idx], dtype=object)
                    if not any(v != 0 for v in xv): continue
                    tot = tot + np.tensordot(xv, Fh, axes=1) * mpq(1, den)
            else:
                gi, idx, mats = blk['data'][0]
                for k, i in enumerate(idx):
                    if x[i] != 0:
                        tot = tot + np.array([[mpq(a_, b_) for (a_, b_) in row] for row in mats[k]], dtype=object) * x[i]
            if any(v != 0 for v in tot.ravel()): bad += 1
    print(f"[data] exact face check: F_k(x_z) = 0 for all {len(blocks)} blocks and {len(zx)} zeros: {bad == 0} ({bad} failures); "
          f"{time.time() - t0:.1f}s", flush=True)
    nbits = max(max(abs(int(v)).bit_length() for gi, idx, Fh, den in blk['data'] for v in Fh.ravel()) for blk in blocks if blk['kind'] == 'COV')
    with open('facedata4m.pkl', 'wb') as fh: pickle.dump(dict(blocks=blocks, primes=primes), fh)
    print(f"saved facedata4m.pkl ({len(blocks)} blocks, COV numerators up to {nbits} bits); {time.time() - t0:.1f}s", flush=True)

if __name__ == "__main__":
    main()
