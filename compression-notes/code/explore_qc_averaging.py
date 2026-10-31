"""Exploration (not in the notes yet): C-dependent averages for the compression conjecture QC.
For copy s let rho_s = local support operator of C at s (reduced operator of C C^* + C^* C on copy s).  Draw w from the
angular central Gaussian with shape rho_s^p (w = rho_s^{p/2} g / |rho_s^{p/2} g|, g standard complex Gaussian) and
estimate E[Delta_s(w)] by Monte Carlo.  QC follows if, for some fixed p, max_s E[Delta_s(w)] >= 0 for every C.
p = 0 is the C-independent average (fails at t = 1/2); p -> -infinity is the smallest-eigenvector rule (fails on a
stored operator).  usage: python explore_qc_averaging.py samples"""
import sys, glob
import numpy as np
from werner import d, q, applyB, contract
sys.path.insert(0, '.')
exec(open('check_qc.py', encoding='utf-8').read().split('def quad_part')[0].split('rng = np.random.default_rng(0)')[1])
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
t = 0.5
rng2 = np.random.default_rng(7)
def local_support(C, n, s):
    M = C @ C.conj().T + C.conj().T @ C
    T = M.reshape([d] * (2 * n))
    keep = s
    letters = 'abcdefghijklmnopqrstuvwxyz'
    inp = [letters[i] for i in range(2 * n)]
    for k in range(n):
        if k != keep: inp[n + k] = inp[k]
    out = inp[keep] + inp[n + keep]
    return np.einsum(''.join(inp) + '->' + out, T)
def avg_delta(C, n, s, p, X):
    rs = local_support(C, n, s); ev, U = np.linalg.eigh(rs); ev = np.maximum(ev, 1e-14 * ev.max())
    S = U @ np.diag(ev ** (p / 2)) @ U.conj().T
    vals = []
    for _ in range(NS):
        g = rng2.standard_normal(d) + 1j * rng2.standard_normal(d); w = S @ g
        vals.append(delta(C, n, s, w, t, X))
    vals = np.array(vals); return vals.mean(), vals.std() / np.sqrt(NS)
def randC(n):
    N = d ** n; Y = rng2.standard_normal((N, 2)) + 1j * rng2.standard_normal((N, 2)); Z = rng2.standard_normal((N, 2)) + 1j * rng2.standard_normal((N, 2))
    return Y @ Z.conj().T
cases = []
for f in sorted(glob.glob('../data/rule_*_n2.npy')) + sorted(glob.glob('../data/qc_adversarial_n3.npy')):
    C = np.load(f); n = int(round(np.log(C.shape[0]) / np.log(d))); cases.append((f.split('/')[-1], C, n))
for k in range(3): cases.append((f"random n=2 #{k}", randC(2), 2))
for k in range(2): cases.append((f"random n=3 #{k}", randC(3), 3))
ps = [0.0, -0.5, -1.0, -2.0, -4.0]
for name, C, n in cases:
    C = C / np.linalg.norm(C); qn = q(C, n, t)
    Xs = [blocks(C, n, s) for s in range(n)]
    row = []
    for p in ps:
        best = max(avg_delta(C, n, s, p, Xs[s])[0] for s in range(n))
        row.append(best)
    print(f"{name}: q_n = {qn:+.4f}; max_s E[Delta] for p = {ps}: " + ", ".join(f"{v:+.4f}" for v in row), flush=True)
