"""The families of the type I forms with a letter 1 at one marked copy, and the forms they must cover.

Sizes (q0, t, i, tau'); each is fixed at a value below its symbolic minimum (tfamlib.SMIN: q0 >= 4, t >= 4, i >= 2,
tau' >= 2) or symbolic.  Fixed values: q0 in -1..3; t in -1..3 (case MT) or -2..3 (case MU); i in 0, 1; tau' in 0, 1
(MT sym), 1 (MU sym), none (anti, where tau' >= 2 always).  A form (k, j, tau, case, mode, i) is *needed* when
k >= 30, 3 <= j <= k - 2, 1 <= tau <= j and its twisted space is nonzero; in the sizes, for case MT: j = q0 + tau' +
1 + i >= 3, t + i >= 2, q0 + i >= 0; for case MU: j = q0 + tau' + i >= 3, t + i >= 1, q0 + i >= 0, tau' >= 1.
A family is listed when it contains a needed form; the conditions are monotone in the sizes, so it suffices to test
the member with every symbolic size at its minimum + 100 (k >= 30 there).
python famlist.py CASE MODE NSYM_MIN NSYM_MAX  ->  one family per line, as accepted by tprove.py"""
import sys
import tfamlib as T
import tlib

QF = (-1, 0, 1, 2, 3)
TF = {'MT': (-1, 0, 1, 2, 3), 'MU': (-2, -1, 0, 1, 2, 3)}
IF = (0, 1)
PF = {('MT', 'sym'): (0, 1), ('MU', 'sym'): (1,), ('MT', 'anti'): (), ('MU', 'anti'): ()}


def needed_sizes(case, q0, t, i, tp):
    if case == 'MT':
        return q0 + tp + 1 + i >= 3 and t + i >= 2 and q0 + i >= 0 and tp >= 0
    return q0 + tp + i >= 3 and t + i >= 1 and q0 + i >= 0 and tp >= 1


def nonzero(case, mode, q0, t, i, tp):
    return bool(T.present_at(case, mode, {'q0': q0, 't': t, 'i': i, 'tp': tp}))


def families(case, mode, nmin=1, nmax=4):
    out = []
    for q0 in QF + (None,):
        for t in TF[case] + (None,):
            for i in IF + (None,):
                for tp in PF[(case, mode)] + (None,):
                    fam = (q0, t, i, tp)
                    ns = sum(v is None for v in fam)
                    if not (nmin <= ns <= nmax):
                        continue
                    far = [T.SMIN[nm] + 100 if v is None else v for v, nm in zip(fam, T.SIZES)]
                    if not needed_sizes(case, *far):
                        continue
                    if not nonzero(case, mode, *far):
                        continue
                    out.append((ns, fam))
    # cheapest first: fewer symbolic sizes; then i fixed before i symbolic, tau' fixed before symbolic
    out.sort(key=lambda z: (z[0], z[1][2] is None, z[1][3] is None, z[1][0] is None, z[1][1] is None))
    return [fam for _, fam in out]


def spec(case, mode, fam):
    return f'{case}/{mode}/' + ','.join('n' if v is None else str(v) for v in fam)


if __name__ == '__main__':
    case, mode = sys.argv[1], sys.argv[2]
    nmin, nmax = int(sys.argv[3]), int(sys.argv[4])
    for fam in families(case, mode, nmin, nmax):
        print(spec(case, mode, fam))
