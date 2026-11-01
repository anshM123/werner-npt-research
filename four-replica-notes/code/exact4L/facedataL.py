"""Exact cross data of the lifted problem at n = 4 (stage 2 of the lifted constraint data).

Faces (read only): N+_k from ../exact4/face4.pkl (t = +1 zeros), N-_k from ../exact4m/face4m.pkl (t = -1 zeros); both use
the same block list (X groups, COV class blocks with identical compressions U, checked here).
COV block k = (xy, mu, sigma), N = [N+ | N-] (2m x d, d = d+ + d-, rows split as N_r = rows r*m..(r+1)*m): for every
parameter i of group g
    Hhat_{k,i} = H0 + H0^T,   H0 = sum_{r,s} (Gt U N_r)^T Ytil_rs,i (U N_s),      H_k(e_i) = Hhat_{k,i} / (2 dT^n dX_g),
i.e. H = N^T W N with W = W0 + W0^T as in ../exact4/facedata4.py, so
    Hhat = [[F+hat, Ghat], [Ghat^T, F-hat]],   Ghat = N+^T W0 N- + (N-^T W0 N+)^T   (cross data G_k = N+^T W_k N-).
H0 is computed modulo 27-bit primes from the float64-exact tree images and recovered by CRT with the a-priori bound
    |Hhat| <= 2 sum_{r,s} max|Ytil_rs| colsum(Gt U N_r) colsum(U N_s),   max|Ytil| <= |orbit| tau^n Xmax_g.
The work is done in three passes (xy = uu, uv, vv) to bound the memory.
X block of group g: H_g(e_k) = N^T Y_k N exactly (flint), N = [N+_g | N-_g].
Checks (all exact):
  * the F+ and F- parts equal the certified face data ../exact4/facedata4.pkl and ../exact4m/facedata4m.pkl entry by
    entry (blocks / groups with identically zero data are absent there and must vanish here);
  * G_k(x_z) = 0 for all blocks at all 27 zeros (W(x_z) N+ = 0 at the 11 t = +1 zeros, W(x_z) N- = 0 at the 16 t = -1
    zeros, so N+^T W(x_z) N- = 0 at both kinds); F+_k(x_z) = 0 at the t = +1 zeros, F-_k(x_z) = 0 at the t = -1 zeros;
    the numbers of blocks with F+_k(x_z) != 0 at t = -1 zeros (and F-_k(x_z) != 0 at t = +1 zeros) are reported (they
    need not vanish).
Saves facedataL.pkl: blocks = list of dict(name, kind, dp, dm, data = [(gi, idx, Ghat (K, dp, dm) object ints, den)]),
one entry per block with dp + dm > 0 (data empty when G vanishes identically); G_k(e_i) = Ghat[kk] / den.
usage: python facedataL.py"""
import sys, os, time, pickle, math, itertools, numpy as np, flint, gmpy2
from gmpy2 import mpq
from core4 import *
from engine4 import *

def _primes_below(B, k):
    out = []; q = B
    while len(out) < k:
        q -= 1
        if gmpy2.is_prime(q): out.append(q)
    return out
PRIMES = _primes_below(2 ** 27, 16)        # same primes as ../exact4/points4.py

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

def eq_arrays(a, b):
    return a.shape == b.shape and all(u == v for u, v in zip(a.ravel(), b.ravel()))

def frac_mats_to_int(mats):
    """list (over parameters) of r x r lists of (num, den) -> (object int array (K, r, r), common denominator)"""
    d = 1
    for M in mats:
        for row in M:
            for (a, q) in row: d = d * q // math.gcd(d, q)
    arr = np.array([[[a * (d // q) for (a, q) in row] for row in M] for M in mats], dtype=object)
    return arr, d

def main():
    t0 = time.time(); n = 4
    A, B1, A3 = site_algebras()
    vs = VarSpace(A, n, log=lambda s: None); N = vs.N
    fp = pickle.load(open('../exact4/face4.pkl', 'rb')); fm = pickle.load(open('../exact4m/face4m.pkl', 'rb'))
    assert fp['LX'] == fm['LX'] and fp['dT'] == fm['dT']
    for k in fp['U']:
        a, b = fp['U'][k], fm['U'][k]
        assert len(a) == len(b) and all(la == lb and da == db and eq_arrays(Ua, Ub) for (la, da, Ua), (lb, db, Ub) in zip(a, b)), \
            f"class compressions differ for {k}"
    assert set(fp['cker']) == set(fm['cker']) and set(fp['xker']) == set(fm['xker']) == set(vs.order)
    print(f"[dataL] faces: t=+1 zeros {fp['names']}; t=-1 zeros {fm['names']}; class compressions U identical in both; "
          f"{time.time() - t0:.1f}s", flush=True)
    ref = {}
    for tag, fn in (('+', '../exact4/facedata4.pkl'), ('-', '../exact4m/facedata4m.pkl')):
        fd = pickle.load(open(fn, 'rb'))
        assert list(fd['primes']) == PRIMES[:len(fd['primes'])]
        dd = {}
        for b in fd['blocks']:
            if b['kind'] == 'COV': dd[b['name']] = {gi: (Fh, den) for gi, idx, Fh, den in b['data']}
            else:
                gi, idx, mats = b['data'][0]; dd[b['name']] = {gi: frac_mats_to_int(mats)}
        ref[tag] = dd; del fd
    nref = {tag: sum(len(v) for v in ref[tag].values()) for tag in ref}
    print(f"[dataL] reference face data loaded: {len(ref['+'])} t=+1 blocks ({nref['+']} block-groups), {len(ref['-'])} t=-1 "
          f"blocks ({nref['-']} block-groups); {time.time() - t0:.1f}s", flush=True)
    matched = {'+': 0, '-': 0}; bad = {'+': 0, '-': 0}
    mus = sorted_patterns(len(B1.blocks), n)
    Gts = {mu: np.array(fm_to_int(fm_kron([B1.blocks[q]['G'] for q in mu]), 1), dtype=object) for mu in mus}
    maps = {}
    for (xy, r, s_, c, P) in cov_entries():
        Pm = PERM[P]
        Ti, dT = site_map_int(A, B1, lambda X, Pm=Pm: realign1(Pm @ X @ Pm.T))
        maps[(xy, r, s_)] = (c, Ti, dT)
    dT = fp['dT']; assert all(m[2] == dT for m in maps.values())
    tau = max(rowsum_bound(m[1]) for m in maps.values())
    gmax = max(len(vs.groups[k]['orbit']) * vs.groups[k]['Xmax'] for k in vs.order) * tau ** n
    out = []          # saved blocks
    def compare(tag, name, gi, arr, den):
        """exact comparison of a recomputed face part with the reference (absent there <=> identically zero)"""
        nz = any(v != 0 for v in arr.ravel())
        got = ref[tag].get(name, {}).get(gi)
        if not nz:
            if got is not None: bad[tag] += 1
            return
        if got is None: bad[tag] += 1; return
        Fr, dr = got
        # same rational values: arr / den == Fr / dr
        if Fr.shape != arr.shape or not all(u * dr == v * den for u, v in zip(arr.ravel(), Fr.ravel())): bad[tag] += 1
        else: matched[tag] += 1
    # ------------------------------------------------------------------ COV blocks, one pass per xy
    maxbits = 0; nprimes_used = set(); stats = []
    for xy in ('uu', 'uv', 'vv'):
        kept = []
        for (xy_, mu), cl in fp['U'].items():
            if xy_ != xy: continue
            for lab, dsig, U in cl:
                Np = fp['cker'][(xy, mu, lab)]; Nm = fm['cker'][(xy, mu, lab)]; m = U.shape[1]
                dp, dm = Np.shape[1], Nm.shape[1]
                if dp + dm == 0: continue
                Nc = np.hstack([Np, Nm]) if dp and dm else (Np if dp else Nm)
                NL = [Nc[:m], Nc[m:]]
                L = [Gts[mu].dot(U).dot(NL[r]) for r in range(2)]
                R = [U.dot(NL[s_]) for s_ in range(2)]
                kept.append(dict(xy=xy, mu=mu, lab=lab, m=m, dp=dp, dm=dm, r=dp + dm, L=L, R=R))
        kidx = {}
        for bi, b in enumerate(kept): kidx.setdefault(b['mu'], []).append(bi)
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
        nprimes_used.add(len(primes))
        print(f"[dataL] pass {xy}: {len(kept)} COV blocks (max d+ + d- = {max(b['r'] for b in kept)}); |Hhat| <= 2^{fb.bit_length()}: "
              f"{len(primes)} primes; {time.time() - t0:.1f}s", flush=True)
        def modp(M, p): return np.array([[int(v) % p for v in row] for row in M], dtype=np.int64).reshape(M.shape)
        Lst = {}; Rst = {}; offs = {}
        for mu, bis in kidx.items():
            o = [0]
            for bi in bis: o.append(o[-1] + kept[bi]['r'])
            offs[mu] = o
            for p in primes:
                for r in range(2):
                    Lst[(mu, r, p)] = modp(np.hstack([kept[bi]['L'][r] for bi in bis]), p)
                    for bi in bis: Rst[(bi, r, p)] = modp(kept[bi]['R'][r], p)
        acc = {}
        for r in range(2):
            for s_ in range(2):
                c, Ti, _ = maps[(xy, r, s_)]
                for key in vs.order:
                    if key[0] != c: continue
                    g = vs.groups[key]
                    level = tree(vs.pair_float(key), key[1], Ti, len(B1.blocks))
                    bnd = g['Xmax'] * tau ** n * len(g['orbit'])
                    for mu in mus:
                        if mu not in kidx: continue
                        Y = assemble(level, g['orbit'], mu)
                        if Y is None: continue
                        check_exact(Y, bnd)
                        gs = tuple(B1.blocks[q]['f'] for q in mu)
                        Yi = pair_to_mat(Y, gs).astype(np.int64)
                        o = offs[mu]
                        for p in primes:
                            Yp = Yi % p
                            T1 = np.matmul(Lst[(mu, r, p)].T, Yp) % p          # (K, sum d, Dm); sums of <= 81 terms < 2^54
                            for j, bi in enumerate(kidx[mu]):
                                T2 = np.matmul(T1[:, o[j]:o[j + 1], :], Rst[(bi, s_, p)]) % p
                                kk = (bi, key)
                                if kk not in acc:
                                    acc[kk] = {q: np.zeros((g['K'], kept[bi]['r'], kept[bi]['r']), dtype=np.int64) for q in primes}
                                acc[kk][p] = (acc[kk][p] + T2) % p
                print(f"   entry {xy}{'uv'[r]}{'uv'[s_]} (family {c}); {time.time() - t0:.1f}s", flush=True)
        # symmetrise, CRT, split into F+, G, F-
        for bi, b in enumerate(kept):
            name = f"COV{b['xy']}{b['mu']}{b['lab']}"; dp, dm = b['dp'], b['dm']
            data = []
            for key in vs.order:
                kk = (bi, key)
                if kk not in acc: continue
                res = [(acc[kk][p] + acc[kk][p].transpose(0, 2, 1)) % p for p in primes]
                Hh = crt_array(res, primes)
                del acc[kk]
                g = vs.groups[key]; gi = vs.order.index(key); den = 2 * dT ** n * g['dX']
                if any(v != 0 for v in Hh.ravel()):
                    maxbits = max(maxbits, max(abs(int(v)).bit_length() for v in Hh.ravel()))
                if dp: compare('+', name, gi, Hh[:, :dp, :dp], den)
                if dm: compare('-', name, gi, Hh[:, dp:, dp:], den)
                if dp and dm:
                    Gh = Hh[:, :dp, dp:]
                    if any(v != 0 for v in Gh.ravel()): data.append((gi, list(g['idx']), np.array(Gh, dtype=object), den))
            # a reference block-group never produced here is caught by the final count (matched == number of references)
            out.append(dict(name=name, kind='COV', dp=dp, dm=dm, data=data))
        acc.clear(); del acc, Lst, Rst
        print(f"[dataL] pass {xy} done: running exact comparison with the saved face data: t=+1 {matched['+']} block-groups "
              f"equal, {bad['+']} mismatches; t=-1 {matched['-']} equal, {bad['-']} mismatches; {time.time() - t0:.1f}s", flush=True)
    # ------------------------------------------------------------------ X blocks
    for key in vs.order:
        Np, Nm = fp['xker'][key], fm['xker'][key]; dp, dm = Np.shape[1], Nm.shape[1]
        if dp + dm == 0: continue
        g = vs.groups[key]; t = g['t']; gi = vs.order.index(key); name = f"X{key[0]}{key[1]}"
        Nc = np.hstack([Np, Nm]) if dp and dm else (Np if dp else Nm)
        Nq = flint.fmpq_mat(t, dp + dm, [int(v) for v in Nc.ravel()])
        Hs = [Nq.transpose() * g['Ys'][k] * Nq for k in range(g['K'])]
        den = 1
        for H in Hs:
            for a in range(dp + dm):
                for b_ in range(dp + dm): den = den * int(H[a, b_].q) // math.gcd(den, int(H[a, b_].q))
        Hh = np.array([[[int((H[a, b_] * den).p) for b_ in range(dp + dm)] for a in range(dp + dm)] for H in Hs], dtype=object)
        assert all(Hh[k, a, b_] * 1 == Hh[k, b_, a] for k in range(g['K']) for a in range(dp + dm) for b_ in range(dp + dm))
        if dp: compare('+', name, gi, Hh[:, :dp, :dp], den)
        if dm: compare('-', name, gi, Hh[:, dp:, dp:], den)
        data = []
        if dp and dm:
            Gh = Hh[:, :dp, dp:]
            if any(v != 0 for v in Gh.ravel()): data.append((gi, list(g['idx']), np.array(Gh, dtype=object), den))
        out.append(dict(name=name, kind='X', dp=dp, dm=dm, data=data))
    allok = (bad['+'] == 0 and bad['-'] == 0 and matched['+'] == nref['+'] and matched['-'] == nref['-'])
    print(f"[dataL] exact comparison of the F+ / F- parts with ../exact4/facedata4.pkl and ../exact4m/facedata4m.pkl: "
          f"t=+1 {matched['+']} of {nref['+']} block-groups equal ({bad['+']} mismatches), t=-1 {matched['-']} of {nref['-']} "
          f"equal ({bad['-']} mismatches): {allok}; COV numerators up to {maxbits} bits; {time.time() - t0:.1f}s", flush=True)
    nG = sum(1 for b in out if b['data']); nGent = sum(Gh.size for b in out for (_, _, Gh, _) in b['data'])
    print(f"[dataL] cross data: {len(out)} blocks with d+ + d- > 0 ({sum(1 for b in out if b['kind'] == 'COV')} COV, "
          f"{sum(1 for b in out if b['kind'] == 'X')} X); {nG} with G not identically zero; {nGent} numerator entries; "
          f"{time.time() - t0:.1f}s", flush=True)
    # ------------------------------------------------------------------ exact checks at the zeros
    def load_zero(path, t):
        """zero as an integer vector a with x = a / lambda (lambda > 0 the lcm of the denominators)"""
        rec = pickle.load(open(path, 'rb'))
        if t == 1: assert rec['in_space'] and rec['eq_ok'] and rec['mom_ok'] and rec['obj'] == (0, 1)
        else: assert rec.get('t') == -1 and rec['kind'] != 'ent' and rec['q4'] == (0, 1) and rec['zero_ok']
        lam = 1
        for (p_, q_) in rec['x']: lam = lam * q_ // math.gcd(lam, q_)
        return [p_ * (lam // q_) for (p_, q_) in rec['x']]
    zp = [(nm, load_zero(f'../exact4/pts/{nm}.pkl', 1)) for nm in fp['names']]
    zm = [(nm, load_zero(f'../exact4m/pts/{nm}.pkl', -1)) for nm in fm['names']]
    def combo(items, a):
        """lambda * L * sum_g (1/den_g) sum_k x_k M_k for items = [(idx, M (K, ., .), den_g)], L = lcm of the den_g
        (integer arithmetic; zero iff the rational block value is zero)"""
        L = 1
        for idx, M, den in items: L = L * den // math.gcd(L, den)
        tot = None
        for idx, M, den in items:
            av = np.array([a[i] for i in idx], dtype=object)
            if not any(v != 0 for v in av): continue
            T = np.tensordot(av, M, axes=1) * (L // den)
            tot = T if tot is None else tot + T
        return tot
    def evalG(b, a): return combo([(idx, Gh, den) for gi, idx, Gh, den in b['data']], a)
    def evalref(tag, name, a):
        return combo([(vs.groups[vs.order[gi]]['idx'], Fr, dr) for gi, (Fr, dr) in ref[tag].get(name, {}).items()], a)
    nonzero = lambda T: T is not None and any(v != 0 for v in T.ravel())
    badG = 0; badFp = 0; badFm = 0; crossp = set(); crossm = set()
    for nm, x in zp + zm:
        for b in out:
            if nonzero(evalG(b, x)): badG += 1
    for nm, x in zp:
        for name in ref['+']:
            if nonzero(evalref('+', name, x)): badFp += 1
        for name in ref['-']:
            if nonzero(evalref('-', name, x)): crossm.add(name)
    for nm, x in zm:
        for name in ref['-']:
            if nonzero(evalref('-', name, x)): badFm += 1
        for name in ref['+']:
            if nonzero(evalref('+', name, x)): crossp.add(name)
    zok = badG == 0 and badFp == 0 and badFm == 0
    print(f"[dataL] exact checks at the zeros: G_k(x_z) = 0 for all {len(out)} blocks at all {len(zp) + len(zm)} zeros "
          f"({badG} failures); F+_k(x_z) = 0 at the {len(zp)} t=+1 zeros ({badFp} failures); F-_k(x_z) = 0 at the {len(zm)} "
          f"t=-1 zeros ({badFm} failures): {zok}.  Not required and not zero: F+_k(x_z) != 0 at some t=-1 zero for "
          f"{len(crossp)} of {len(ref['+'])} t=+1 face blocks, F-_k(x_z) != 0 at some t=+1 zero for {len(crossm)} of "
          f"{len(ref['-'])} t=-1 face blocks; {time.time() - t0:.1f}s", flush=True)
    assert allok and zok
    with open('facedataL.pkl', 'wb') as fh:
        pickle.dump(dict(blocks=out, primes=PRIMES[:max(nprimes_used)], zeros_plus=fp['names'], zeros_minus=fm['names']), fh)
    print(f"saved facedataL.pkl ({len(out)} blocks); {time.time() - t0:.1f}s", flush=True)

if __name__ == "__main__":
    main()
