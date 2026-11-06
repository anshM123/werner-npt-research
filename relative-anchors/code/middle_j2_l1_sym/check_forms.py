"""The certified symbolic forms (saved by batch.py) against the same reduction in exact rational arithmetic at points.

python check_forms.py REGIME [REGIME ...]     (REGIME as in batch.py, e.g. 9+,2)
For each family with a form: evaluate the saved W-block of S' (strings of the numerators and denominators over
Q(q, x...)) at q in {1/9, 4, 49/25} and x in {0, 1, 5} (each symbolic level), and compare with r2v.reduce_family over
the rationals at the integer sizes o + x (the same bases).  The rational-arithmetic reduction itself is compared with the
explicit block from the Lemma B builder by check_points.py (also at k = 11, inside the regimes).  PASS iff all equal."""
import sys, time, pickle, itertools, re
from fractions import Fraction as Fr
import r2v

t0 = time.time()


def ev(s, env):
    """exact evaluation: every integer literal that is not an exponent or part of a name becomes a Fraction"""
    s = re.sub(r'(?<![*\w])(\d+)', r'Fr(\1)', s.replace('^', '**'))
    return eval(s, {'__builtins__': {}, 'Fr': Fr}, env)


def parse(sp):
    return (None, int(sp[:-1])) if sp.endswith('+') else (int(sp), None)


if __name__ == '__main__':
    ok = True
    for reg in sys.argv[1:]:
        s0, s1 = reg.split(',')
        (f0, o0), (f1, o1) = parse(s0), parse(s1)
        forms = pickle.load(open(f'logs/forms_{s0}_{s1}.pkl', 'rb'))
        sample = (o0 + 3 if f0 is None else f0, o1 + 3 if f1 is None else f1)
        xs0 = [0, 1, 5] if f0 is None else [None]
        xs1 = [0, 1, 5] if f1 is None else [None]
        nchk = 0
        for F in r2v.families():
            form = forms.get(F.name)
            if form is None:
                continue
            bases = r2v.choose(F, sample)
            for q in (Fr(1, 9), Fr(4), Fr(49, 25)):
                for x0, x1 in itertools.product(xs0, xs1):
                    u0 = f0 if f0 is not None else o0 + x0
                    u1 = f1 if f1 is not None else o1 + x1
                    env = {'q': q}
                    if x0 is not None:
                        env['x0'] = Fr(x0)
                    if x1 is not None:
                        env['x1'] = Fr(x1)
                    Ssym = [[ev(n, env) / ev(d, env) for (n, d) in row] for row in form['S']]
                    red = r2v.reduce_family(F, r2v.Ctx(Fr(1), q, Fr(u0), Fr(u1)), bases, full=False)
                    eq = Ssym == red['S']
                    nchk += 1
                    if not eq:
                        print(f'  MISMATCH {F.name} at q = {q}, sizes ({u0}, {u1})')
                    ok &= eq
        print(f'({s0}, {s1}): {nchk} comparisons of saved forms with the rational reduction [{time.time() - t0:.1f}s]',
              flush=True)
    print('PASS' if ok else 'FAIL')
