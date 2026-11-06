"""Control for the frame constant of diag_certificate.js (which uses H = 3/4): the same program with H = 4/5 must fail
(the numerical infimum of lambda_min(Ghat_tt)/(1 - theta_t) is 3/4), and with the earlier constant H = 2/3 it passes.
Writes temporary copies next to this file and removes them."""
import os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, 'diag_certificate.js'), encoding='utf-8').read()
assert src.count('sfq(3,4)') >= 4
results = {}
for tag, H in (('4/5', 'sfq(4,5)'), ('2/3', 'sfq(2,3)')):
    path = os.path.join(here, '_control_cert.js')
    open(path, 'w', encoding='utf-8', newline='\n').write(src.replace('sfq(3,4)', H))
    try:
        r = subprocess.run([sys.executable, os.path.join(here, 'run_js.py'), path, 'sfP6Verify'], capture_output=True, text=True)
    finally:
        os.remove(path)
    out = r.stdout + r.stderr
    passed = r.returncode == 0 and '"status": "PASS"' in out[:200]
    failed_as_expected = 'coefficient positivity failure' in out
    results[tag] = (passed, failed_as_expected)
    k = out.find('Error: P6 coefficient')
    print(f'H = {tag}: {"PASS" if passed else "FAIL"}' + (f' ({out[k + 7:out.find(chr(10), k)].strip()})' if k >= 0 else ''))
ok = results['4/5'] == (False, True) and results['2/3'][0]
print('PASS' if ok else 'FAIL')
