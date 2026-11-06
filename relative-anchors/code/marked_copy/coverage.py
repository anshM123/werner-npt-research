"""Coverage of the fixed-k certificates and the certified objects against the Lemma B builder.
(1) For every k = KMIN..KMAX, every j = -1..k and every i = 0..floor((k-1)/2): either the twisted space V_i is empty
    (all twisted vectors have norm 0, by counting) or the certificate run saved a record with status
    'PD on complement', whose kernel vectors are exactly the expected ones: p x h_y0 for each y0 with h_y0 != 0,
    K2 at (j, i) = (k-2, 0), T0 at (k-1, 0), T1 at (k-1, 1).
(2) For k = KMIN..KC: every saved gauge-fixed form, evaluated at sigma in {4, 1/9}, equals the restriction of the
    rescaled block built by j1lib.py (Lemma B) to the kept explicit twisted vectors, entry by entry.
Usage: python coverage.py cert.pkl KMIN KMAX KC"""
import sys, time, pickle
from fractions import Fraction as Fr
import flint
import mlib, explicit

pk, KMIN, KMAX, KC = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
saved = pickle.load(open(pk, 'rb'))
T0 = time.time()
ok = True
ncert = nempty = 0
for k in range(KMIN, KMAX + 1):
    for j in range(-1, k + 1):
        for i in range(0, (k - 1) // 2 + 1):
            C = mlib.Calc(k, j, i, sigma=Fr(1))
            present = [nm for nm in mlib.NAMES if mlib.norm_of(C, nm) != 0]
            if not present:
                nempty += 1
                if (k, j, i) in saved:
                    ok = False
                    print('SAVED BUT EMPTY', k, j, i)
                continue
            rec = saved.get((k, j, i))
            if rec is None or rec.get('status') != 'PD on complement':
                ok = False
                print('NOT CERTIFIED', k, j, i, rec and rec.get('status'))
                continue
            if rec['names'] != present:
                ok = False
                print('NAMES DIFFER', k, j, i)
            exp = []
            for y0 in (0, 1):
                if ('O%d' % y0) in present:
                    exp.append('p x h%d' % y0)
            if j == k - 2 and i == 0:
                exp.append('K2')
            if j == k - 1 and i == 0:
                exp.append('T0')
            if j == k - 1 and i == 1:
                exp.append('T1')
            if rec['kernel'] != exp:
                ok = False
                print('KERNEL LIST DIFFERS', k, j, i, rec['kernel'], exp)
            ncert += 1
print(f'(1) k = {KMIN}..{KMAX}: {ncert} certified reduced forms, {nempty} empty twisted spaces, nothing missing',
      flush=True)


def ev(coeffs, x):
    s = Fr(0)
    for c in reversed(coeffs):
        s = s * x + Fr(c)
    return s


ncmp = 0
for k in range(KMIN, KC + 1):
    for w0 in (Fr(2), Fr(1, 3)):
        sig = w0 * w0
        w = [w0] + [Fr(1)] * (k - 1)
        for j in range(-1, k + 1):
            st, M = explicit.block(k, j, w)
            for i in range(0, (k - 1) // 2 + 1):
                rec = saved.get((k, j, i))
                if rec is None:
                    continue
                keep = [nm for nm in rec['names'] if nm not in rec['gauge']]
                st2, vecs = explicit.twisted(k, j, i)
                for a, na in enumerate(keep):
                    for b, nb in enumerate(keep):
                        num, den = rec['form'][a][b]
                        val = ev(num, sig) / ev(den, sig)
                        ncmp += 1
                        if val != explicit.form(M, vecs[na], vecs[nb]):
                            ok = False
                            print('SAVED FORM DIFFERS FROM LEMMA B', k, j, i, na, nb, sig)
    print(f'(2) k={k}: done ({time.time() - T0:.0f}s)', flush=True)
print(f'(2) {ncmp} entries of saved gauge-fixed forms equal the Lemma B projections')
print('PASS' if ok else 'FAIL', f'{time.time() - T0:.0f}s')
