"""Exact checks of the certified blocks against the two builders.
 (a) Lemma B: for every family and every configuration, the block saved by the certificate (built over Q(q) with
     ylib.py) evaluated at rational squares q0 equals the rescaled block built by j1lib.py (Lemma B) at the weights
     w = sqrt(s) (j1lib gives 2^Z Y; ylib gives Y; Z = k - j).
 (b) Definitions: the block built by ylib.py over Q at a rational point (point ring) equals the saved block there.
 (c) At general weights (distinct squared weights; lists for k <= 7) the two builders agree for every family, also with the B-part
     scaled by 1001/1000 (the builder of the negative controls).
 (d) Homogeneity: scaling every squared weight by c multiplies the rescaled block by c^(j + 1 - |T|).
python check_builders.py FAMILIES [k=5]    (FAMILIES comma separated, e.g. L3,T31)
Prints PASS or FAIL on the last line."""
import sys, time
from fractions import Fraction as Fr
import flint
import ylib
import k5lib as K
from saved import load, block_at, cert_path

t0 = time.time()
args = [a for a in sys.argv[1:] if '=' not in a]
opts = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
k = int(opts.get('k', 5))
fams = args[0].split(',')
ok = True


def log(*a):
    print(*a, flush=True)


def ylib_point(j, T, s, bscale=1):
    sv = [flint.fmpq(Fr(x).numerator, Fr(x).denominator) for x in s]
    st, M = K.build_ylib(k, j, T, sv, flint.fmpq(0), flint.fmpq(1),
                         bscale=bscale if bscale == 1 else flint.fmpq(Fr(bscale).numerator, Fr(bscale).denominator))
    return st, [[Fr(int(x.p), int(x.q)) for x in row] for row in M]


for f in fams:
    j, m = K.fam(f)
    T = tuple(range(m))
    na = nb = 0
    for (a, b) in (K.rep_configs(k, m) if opts.get('reps') == '1' else K.configs(k, m)):
        data = load(cert_path(f, a, b, k))
        for q0 in (Fr(4, 9), Fr(25, 4), Fr(1, 10 ** 4)):
            s = K.config_s(k, m, a, b, q0, Fr(1))
            w = [K.rsqrt(x) for x in s]
            st, Mj = K.build_j1(k, j, T, w)
            assert st == data['states']
            eq = block_at(data, q0) == Mj
            ok &= eq
            na += 1
            if not eq:
                log(f'  (a) MISMATCH {f} {(a, b)} q = {q0}')
        q0 = Fr(9, 16)
        s = K.config_s(k, m, a, b, q0, Fr(1))
        st, My = ylib_point(j, T, s)
        eq = block_at(data, q0) == My
        ok &= eq
        nb += 1
        if not eq:
            log(f'  (b) MISMATCH {f} {(a, b)} q = {q0}')
    log(f'{f}: (a) saved block = Lemma B block at {na} (configuration, q) pairs; (b) = ylib point block at {nb} '
        f'configurations: {ok} [{time.time() - t0:.1f}s]')
    # (c) general weights, both scalings
    for w in ([Fr(1), Fr(2, 3), Fr(3, 4), Fr(5, 6), Fr(7, 2), Fr(9, 7), Fr(11, 5), Fr(13, 4)][:k],
              [Fr(3, 2), Fr(1, 5), Fr(4, 3), Fr(1), Fr(5, 7), Fr(8, 3), Fr(6, 5), Fr(9, 4)][:k]):
        assert len(w) == k, 'general weight vectors are defined for k <= 8'
        s = [x * x for x in w]
        for bs in (Fr(1), Fr(1001, 1000)):
            st, Mj = K.build_j1(k, j, T, w, bscale=bs)
            st2, My = ylib_point(j, T, s, bscale=bs)
            eq = (st == st2) and Mj == My
            ok &= eq
            log(f'  (c) general weights s = {[str(x) for x in s]}, bscale {bs}: ylib = Lemma B: {eq}')
        # (d) homogeneity
        c = Fr(4)
        st, M1 = K.build_j1(k, j, T, w)
        st, M2 = K.build_j1(k, j, T, [2 * x for x in w])
        d = j + 1 - m
        eq = all(M2[i][l] == c ** d * M1[i][l] for i in range(len(st)) for l in range(len(st)))
        ok &= eq
        log(f'  (d) homogeneous of degree {d} in the squared weights: {eq}')
log(f'time {time.time() - t0:.1f}s')
log('PASS' if ok else 'FAIL')
