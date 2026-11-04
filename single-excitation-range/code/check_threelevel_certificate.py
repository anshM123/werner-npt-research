"""Checks for Proposition 7 (single-excitation-range/README.md): the three-level certificate.
(a) the exact one-parameter family: 2Kt equals, for every mu, the pattern + string + child expression (random mu, random
    POVM, n = 4, 5, every sector), via the identity  T - alpha^*alpha = Vac - <Theta, S Theta>;
(b) the heavy-copy family where Proposition 6 fails (n >= 7): on the invariant subspace in which every excited light
    copy carries letter 1 (it contains the minimiser) the three-level certificate is feasible for n = 5..10;
(c) general real POVMs at n = 5, c = 3 (h = 2);
(d) the balanced certificate (mu = 1/2, theta = 0: only child splits of the pattern part) at the same POVMs and at one
    more heavy-copy POVM;
(e) at a near rank-one real frame the splits must be matrices: with k_{S,u} = w_{S,u} R_S (scalar weights, optimised)
    the child conditions fail, while the matrix-valued balanced certificate exists;
every solution is re-checked with numpy eigenvalues, and the guard lambda_min(2Kt) >= (n + 1) * slack is verified."""
import numpy as np, itertools, warnings
from threelevel_certificate import Model, GeneralModel, gauge, certificate
warnings.filterwarnings('ignore')
rng = np.random.default_rng(7)
ok = True
# (a)
worst = 0
for n in (4, 5):
    for c in range(1, n + 1):
        M = GeneralModel(gauge([rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)])); o = M.ops(c)
        worst = max(worst, np.abs((o['T'] - o['AA']) - (o['Vac'] - o['ThS'])).max(), np.abs(sum(o['Pi'].values()) - o['ThS']).max())
print(f"(a) identity T - AA = Vac - <Theta,S Theta> and its split over children: largest error {worst:.1e}")
ok &= worst < 1e-12
# (b)
for N, d in [(4, 0.05), (6, 0.05), (6, 0.2), (7, 0.05), (7, 0.2), (8, 0.2), (9, 0.05)]:
    M = Model(N, d); c = N - 1; lam, mu, st, chk, tmin, tmax = certificate(M, c, three_level=True, real=True)
    l2 = certificate(M, c, three_level=False, real=True)[0]; lm = 2 * np.linalg.eigvalsh(M.ops(c)['K'])[0]
    print(f"(b) heavy family n={N + 1} c={c} d={d}: three-level slack {lam:+.3e} (re-checked {chk:+.3e}, mu {mu:+.3f}); "
          f"Proposition 6 slack {l2:+.3e}; lambda_min(2Kt) - (n+1) slack {lm - (N + 2) * lam:+.3e}", flush=True)
    ok &= chk > -1e-8 and lm - (N + 2) * lam > -1e-7
# (c)
def rz(*s): return rng.normal(size=s)
cases = [("random real", gauge([rz(2, 2) for _ in range(5)]))]
X = [np.outer(rz(2), rz(2)) + 0.05 * rz(2, 2) for _ in range(5)]; cases.append(("near rank-one real frame", gauge(X)))
X = [np.outer(rz(2), rz(2)) + 0.05 * rz(2, 2) for _ in range(5)]; X[0] *= 0.02; X[1] *= 0.02; cases.append(("two nearly unused copies", gauge(X)))
X = [np.eye(2) + 0.1 * rz(2, 2), np.eye(2) + 0.1 * rz(2, 2)] + [0.1 * rz(2, 2) for _ in range(3)]; cases.append(("two heavy copies", gauge(X)))
cases.append(("uniform W-plane", gauge([np.eye(2)] * 5)))
for name, U in cases:
    M = GeneralModel(U); lam, mu, st, chk, tmin, tmax = certificate(M, 3, three_level=True, real=True)
    lm = 2 * np.linalg.eigvalsh(M.ops(3)['K'])[0]
    print(f"(c) n=5 c=3 {name}: three-level slack {lam:+.3e} (re-checked {chk:+.3e}, mu {mu:+.3f}); lambda_min(2Kt) {lm:+.3e}; "
          f"lambda_min(2Kt) - 6 slack {lm - 6 * lam:+.3e}", flush=True)
    ok &= lam > 0 and chk > -1e-8 and lm - 6 * lam > -1e-7
# (d)
X = [np.diag([1.0, 1.0]) + 0.03 * rz(2, 2)] + [0.15 * rz(2, 2) for _ in range(4)]; cases.append(("one heavy copy", gauge(X)))
for name, U in cases:
    M = GeneralModel(U); lam, mu, st, chk, tmin, tmax = certificate(M, 3, three_level=True, real=True, fix_mu=0.5, no_theta=True)
    lm = 2 * np.linalg.eigvalsh(M.ops(3)['K'])[0]
    print(f"(d) n=5 c=3 {name}: balanced certificate (mu = 1/2, theta = 0) slack {lam:+.3e} (re-checked {chk:+.3e}); "
          f"lambda_min(2Kt) - 6 slack {lm - 6 * lam:+.3e}", flush=True)
    ok &= lam > 0 and chk > -1e-8 and lm - 6 * lam > -1e-7
# (e)
import cvxpy as cp
def scalar_split_slack(M, c):
    o = M.ops(c); S_list = o['S']; m = o['m']; n = M.n
    supp = lambda s: frozenset(t for t in range(n) if s[t])
    pats = sorted(set(supp(s) for s in S_list), key=sorted)
    idx = {P: [r * m + j for r in range(2) for j, s in enumerate(S_list) if supp(s) == P] for P in pats}
    R2 = np.real(2 * np.eye(2 * m) + o['Vac'] - o['T'])
    wv = {(P, u): cp.Variable() for P in pats for u in range(n) if u not in P}; lam = cp.Variable(); cons = []
    for P in pats: cons.append(sum(wv[P, u] for u in range(n) if u not in P) == 1)
    for Sp, Pi in o['Pi'].items():
        mem = [u for u in sorted(Sp) if (Sp - {u}) in idx]; ind = sum((idx[Sp - {u}] for u in mem), []); D = len(ind)
        E = 0; off = 0
        for u in mem:
            d = len(idx[Sp - {u}]); Z = np.zeros((D, D)); P = Sp - {u}; Z[off:off + d, off:off + d] = R2[np.ix_(idx[P], idx[P])]
            E = E + wv[P, u] * Z; off += d
        cons.append(np.real(Pi[np.ix_(ind, ind)]) + E - lam * np.eye(D) >> 0)
    cp.Problem(cp.Maximize(lam), cons).solve(solver='CLARABEL'); return lam.value
rng2 = np.random.default_rng(11)
rng2.normal(size=(5, 2, 2))   # skip the "random real" draw of the scan that found this frame
X = [np.outer(rng2.normal(size=2), rng2.normal(size=2)) + 0.05 * rng2.normal(size=(2, 2)) for _ in range(5)]
M = GeneralModel(gauge(X)); ls = scalar_split_slack(M, 3)
lam, mu, st, chk, tmin, tmax = certificate(M, 3, three_level=True, real=True, fix_mu=0.5, no_theta=True)
print(f"(e) n=5 c=3 near rank-one frame: best scalar split slack {ls:+.3e}; matrix-valued balanced certificate slack {lam:+.3e} "
      f"(re-checked {chk:+.3e})", flush=True)
ok &= ls < 0 and lam > 0 and chk > -1e-8
assert ok
print("PASS")
