"""Cross-language check: the energy model of the JavaScript programs (base, norm) against our Python objects (psi_exact):
base == v^T A v and norm == ||R/2||_F^2 exactly, at random rational weights and amplitudes (m = 1).
python crossbase.py SEED N"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'typeI_G_m1_relax'))

import sys, random, json, subprocess, os
import flint
from psi_exact import PsiExact
from test_A30 import parts
from bridge import sphere_point
Q = flint.fmpq
seed, N = int(sys.argv[1]), int(sys.argv[2])
rng = random.Random(seed)
cases = []; ours = []
for it in range(N):
    u = [2, 3, 4, 5][it % 4]
    w = sphere_point(rng, u)
    wt = Q(rng.randint(1, 40), rng.randint(1, 40))
    ps = PsiExact(wt, w)
    Amat, _ = parts(ps)
    n = Amat.nrows(); na = u - 1
    vec = [Q(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(n)]
    alpha = vec[:na]; c = vec[na:na + u]; rho = vec[na + u]; zeta = vec[na + u + 1:]
    a = [sum((ps.Bm[i, r] * alpha[r] for r in range(na)), Q(0)) for i in range(u)]
    base = sum((vec[p] * Amat[p, q] * vec[q] for p in range(n) for q in range(n)), Q(0))
    R = ps.Pi * (ps.Gcal(a, c) * ps.sig + ps.Hfull(rho, zeta)) * ps.Pi
    norm = sum(((R[i, j] / 2) ** 2 for i in range(u) for j in range(u)), Q(0))
    x = [wi * wi for wi in w]
    cases.append({'wt': [str(wt)], 'x': [str(v) for v in x], 'mu': [str(v) for v in x],
                  'A': [str(a[i] / w[i]) for i in range(u)], 'C': [str(c[i] / w[i]) for i in range(u)],
                  'Z': [[str(zeta[i] / w[i]) for i in range(u)]], 'rho': [str(rho)]})
    ours.append((str(base), str(norm)))
json.dump(cases, open('cross_cases.json', 'w'))
js = open('relax_comparison.js', encoding='utf-8').read()
harness = js + r'''
;const __fs=require('fs');const __cases=JSON.parse(__fs.readFileSync('cross_cases.json','utf8'));
const __q=s=>{const [a,b]=s.split('/');return sfq(BigInt(a),BigInt(b??'1'));};
const __out=__cases.map(c=>{const r=sfQMeasureModelGeneralExpanded(c.wt.map(__q),c.x.map(__q),c.mu.map(__q),c.A.map(__q),c.C.map(__q),c.Z.map(z=>z.map(__q)),c.rho.map(__q));return [r.base.str(),r.norm.str()];});
console.log(JSON.stringify(__out));
'''
open('cross_harness.js', 'w', encoding='utf-8').write(harness)
res = subprocess.run(['node', '--stack-size=65500', 'cross_harness.js'], capture_output=True, text=True, timeout=1200)
theirs = json.loads(res.stdout.strip().split('\n')[-1])
nb = sum(1 for o, t in zip(ours, theirs) if Q(o[0]) == Q(t[0]))
nn = sum(1 for o, t in zip(ours, theirs) if Q(o[1]) == Q(t[1]))
print(f'{N} cases: base equal {nb}/{N}; norm equal {nn}/{N}')
print('PASS' if nb == nn == N else 'FAIL')
