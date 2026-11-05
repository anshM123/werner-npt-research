"""Generic sym letter-2 family with tau FIXED: [q0 >= 2, t >= 4, i >= 2, tau = tau0] for tau0 = 1..T0
(each a 3-parameter family over QQ(q0, t, i); same steps as prove_letter2_split2.py plus window and sanity checks).
Runs in one process.  Usage: python prove_tau_loop.py T0 [time_budget_s]"""
import sys, io, time, runpy, contextlib
T0 = int(sys.argv[1]); budget = float(sys.argv[2]) if len(sys.argv) > 2 else 1e9
t_start = time.time()
done = []
for tau0 in range(1, T0 + 1):
    if time.time() - t_start > budget:
        print(f'time budget reached before tau = {tau0}', flush=True)
        break
    t0 = time.time()
    sys.argv = ['prove_letter2_split_tau.py', 'S', 'S', 'S', str(tau0)]
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
    print(f'tau = {tau0}: {"PROVED" if ok else "NOT PROVED"} ({time.time()-t0:.0f}s) {last}', flush=True)
    if not ok:
        print(out[-2000:], flush=True)
        break
    done.append(tau0)
print('proved for tau =', done)
