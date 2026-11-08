# Dual of the Wu-Zou certificate LP (wz_certificate_lp.py) at odd n: find pseudo-profiles (c, P) obeying every linearised
# WZ constraint c_{K+A} + c_{K+B} - c_K <= P_0 - P_K, shadow positivity EXCL(c), EXCL(P) >= 0, c_0 = P_0 = 1, c_[n] = 0,
# and optional valid extras (dim, cs, lmin, uhl), maximising B(c) - A(P).  A positive optimum: no certificate exists.
# that maximise B(c) - A(P).  A positive optimum = the reason no certificate exists.
import itertools, sys, numpy as np
from scipy.optimize import linprog
n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
subs = [frozenset(S) for k in range(n + 1) for S in itertools.combinations(range(n), k)]
idx = {S: i for i, S in enumerate(subs)}; m = len(subs); full = frozenset(range(n))
def incl(R):
    v = np.zeros(m)
    for K in subs:
        if K <= R: v[idx[K]] += (-1) ** len(K) / 2 ** len(R)
    return v
INCL = {R: incl(R) for R in subs}
def excl(S):
    v = np.zeros(m)
    for R in subs:
        if R >= S: v += (-1) ** len(R - S) * INCL[R]
    return v
EXCL = {S: excl(S) for S in subs}
def level(r): return sum(EXCL[S] for S in subs if len(S) == r)
q9 = lambda k: (9 ** k - 1) // 8
B = sum(q9(k) * (3 * level(2 * k + 1) - level(2 * k)) for k in range(1, (n - 1) // 2 + 1))
A = sum(q9(k) * level(2 * k) for k in range(1, n // 2 + 1))
# variables: c (m), P (m)   [P will be constrained symmetric under complement]
nv = 2 * m
Aub, bub = [], []
for K in subs:
    rest = sorted(full - K)
    if len(rest) < 2: continue
    seen = set()
    for r in range(1, len(rest)):
        for Ac in itertools.combinations(rest, r):
            As = frozenset(Ac); Bs = frozenset(rest) - As; key = frozenset([As, Bs])
            if key in seen: continue
            seen.add(key)
            row = np.zeros(nv)
            row[idx[K | As]] += 1; row[idx[K | Bs]] += 1; row[idx[K]] -= 1
            row[m + idx[frozenset()]] -= 1; row[m + idx[K]] += 1      # <= P_0 - P_K
            Aub.append(row); bub.append(0)
for R in subs:
    row = np.zeros(nv); row[:m] = -EXCL[R]; Aub.append(row); bub.append(0)
    row = np.zeros(nv); row[m:] = -EXCL[R]; Aub.append(row); bub.append(0)
Aeq, beq = [], []
row = np.zeros(nv); row[idx[frozenset()]] = 1; Aeq.append(row); beq.append(1)
row = np.zeros(nv); row[m + idx[frozenset()]] = 1; Aeq.append(row); beq.append(1)
row = np.zeros(nv); row[idx[full]] = 1; Aeq.append(row); beq.append(0)
for S in subs:
    T = full - S
    if idx[S] < idx[T]:
        row = np.zeros(nv); row[m + idx[S]] = 1; row[m + idx[T]] = -1; Aeq.append(row); beq.append(0)
EXTRA = sys.argv[2] if len(sys.argv) > 2 else ''
if 'dim' in EXTRA:      # pure-state purity lower bound P_S >= 3^{-min(|S|, n-|S|)}
    for S in subs:
        row = np.zeros(nv); row[m + idx[S]] = -1; Aub.append(row); bub.append(-3.0 ** (-min(len(S), n - len(S))))
if 'cs' in EXTRA:       # cross purity <= purity (Cauchy-Schwarz + AM-GM symmetrisation)
    for S in subs:
        row = np.zeros(nv); row[idx[S]] = 1; row[m + idx[S]] = -1; Aub.append(row); bub.append(0)
if 'lmin' in EXTRA:     # single site: c_s = Tr(rho^u_s rho^v_s) >= lambda_min(rho_s) >= g(P_s) = (1 - sqrt(6P-2))/3 (convex);
    # symmetrised: c_s >= tangent lines of g at points p0 in (1/3, 1/2)
    for p0 in np.linspace(0.334, 0.499, 25):
        g0 = (1 - np.sqrt(6 * p0 - 2)) / 3; dg = -1 / np.sqrt(6 * p0 - 2)    # g'(p) = -(1/3)*(6/(2 sqrt)) = -1/sqrt(6p-2)
        for S in subs:
            if len(S) != 1: continue
            # c_S >= g0 + dg (P_S - p0)   <=>   -c_S + dg P_S <= -g0 + dg p0
            row = np.zeros(nv); row[idx[S]] = -1; row[m + idx[S]] = dg; Aub.append(row); bub.append(-g0 + dg * p0)
if 'uhl' in EXTRA:      # Uhlmann: c_{S^c} = ||X_S||_2^2 <= ||X_S||_1^2 <= F(rho^u_S, rho^v_S)^2 <= rank * Tr(rho^u_S rho^v_S) = r_S c_S
    for S in subs:
        if len(S) == 0 or len(S) == n: continue
        r = 3 ** min(len(S), n - len(S))
        row = np.zeros(nv); row[idx[full - S]] = 1; row[idx[S]] = -r; Aub.append(row); bub.append(0)
obj = np.zeros(nv); obj[:m] = -B; obj[m:] = A        # minimise A(P) - B(c)
res = linprog(obj, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
              bounds=[(0, 1)] * nv, method='highs')
print(f"n={n}: status {res.status}, max B(c) - A(P) = {-res.fun:.6f}")
c = res.x[:m]; P = res.x[m:]
for k in range(n + 1):
    cs = [round(c[idx[S]], 4) for S in subs if len(S) == k]
    Ps = [round(P[idx[S]], 4) for S in subs if len(S) == k]
    print(f"  |S|={k}: c in [{min(cs)}, {max(cs)}] distinct {sorted(set(cs))[:6]};  P distinct {sorted(set(Ps))[:6]}")
np.savez(f'wz_dual_n{n}.npz', c=c, P=P)
