"""Generic sym letter-2 family with ONE parameter fixed (t = value, or q0 = value), the other three symbolic, for a
range of values (3-parameter families; window + sanity + split elimination). Usage: python prove_fixed_loop.py t|q0 lo hi [budget_s]
Runs in one process.  Usage: python prove_tau_loop.py T0 [time_budget_s]"""
import sys, io, time, runpy, contextlib
WHICH = sys.argv[1]; LO, HI = int(sys.argv[2]), int(sys.argv[3]); budget = float(sys.argv[4]) if len(sys.argv) > 4 else 1e9
t_start = time.time()
done = []
for tau0 in range(LO, HI + 1):
    if time.time() - t_start > budget:
        print(f'time budget reached before value = {tau0}', flush=True)
        break
    t0 = time.time()
    sys.argv = ['prove_letter2_split_tau.py'] + (['S', str(tau0), 'S', 'S'] if WHICH == 't' else [str(tau0), 'S', 'S', 'S'])
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            runpy.run_path('prove_letter2_split_tau.py', run_name='__main__')
        out = buf.getvalue()
        ok = 'GENERIC SYM FAMILY PROVED' in out and 'window exact and sanity at sample points: PASS' in out
    except Exception as e:
        out = buf.getvalue() + f'\nEXCEPTION {e!r}'
        ok = False
    last = [l for l in out.splitlines() if 'pivot' in l][-1:]
    print(f'{WHICH} = {tau0}: {"PROVED" if ok else "NOT PROVED"} ({time.time()-t0:.0f}s) {last}', flush=True)
    if not ok:
        print(out[-2000:], flush=True)
        break
    done.append(tau0)
print('proved for', WHICH, '=', done)
