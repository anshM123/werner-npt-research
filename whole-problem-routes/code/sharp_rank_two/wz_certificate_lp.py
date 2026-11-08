"""Can Wu-Zou's inequalities, applied to all traced sets K and all bipartitions A|B of the rest, prove the coupled sharp
bound B_n <= sqrt(A_n(u) A_n(v)) at odd n?  For orthonormal u, v on n copies:
  c_S = Tr(rho^u_S rho^v_S) (= ||Tr_S |u><v| ||_2^2), P_S(w) = Tr(rho^w_S)^2 (= P_{S^c} for pure w), c_0 = P_0 = 1.
WZ (Thm 1.1 + Lemma 5.1):  c_{K+A} + c_{K+B} - c_K <= sqrt(e_K(u) e_K(v)),  e_K = 1 - P_K.
Certificate:  (i)  sum w_{K,AB} (c_{K+A} + c_{K+B} - c_K) - B_n(c) = sum_R l_R excl_uv(R)(c) + kappa c_0   (identity in c),
              (ii) A_n(P) - sum_K W_K (1 - P_K) = sum_S m_S excl_ww(S)(P) + kappa' P_0  (identity in P mod P_S = P_{S^c}),
with W_K = sum_AB w_{K,AB}, all multipliers >= 0.  Then B_n <= sum w sqrt(e e') <= sqrt(sum W e(u) sum W e(v)) <= sqrt(A A')."""
import itertools, sys
import numpy as np
from scipy.optimize import linprog
EXTRA = sys.argv[1] if len(sys.argv) > 1 else 'none'
def run(n):
    subs = [frozenset(S) for k in range(n + 1) for S in itertools.combinations(range(n), k)]
    idx = {S: i for i, S in enumerate(subs)}; m = len(subs); full = frozenset(range(n))
    def incl(R):     # inclusive antisymmetric weight on R as a linear form in c (or P)
        v = np.zeros(m)
        for K in subs:
            if K <= R:
                v[idx[K]] += (-1) ** len(K) / 2 ** len(R)
        return v
    INCL = {R: incl(R) for R in subs}
    def excl(S):
        v = np.zeros(m)
        for R in subs:
            if R >= S:
                v += (-1) ** len(R - S) * INCL[R]
        return v
    EXCL = {S: excl(S) for S in subs}
    def level(r):
        return sum(EXCL[S] for S in subs if len(S) == r)
    q9 = lambda k: (9 ** k - 1) // 8
    B = sum(q9(k) * (3 * level(2 * k + 1) - level(2 * k)) for k in range(1, (n - 1) // 2 + 1))
    A = sum(q9(k) * level(2 * k) for k in range(1, n // 2 + 1))
    # WZ generators
    gens = []
    for K in subs:
        rest = sorted(full - K)
        if len(rest) < 2:
            continue
        seen = set()
        for r in range(1, len(rest)):
            for Acomb in itertools.combinations(rest, r):
                Aset = frozenset(Acomb); Bset = frozenset(rest) - Aset
                key = frozenset([Aset, Bset])
                if key in seen:
                    continue
                seen.add(key)
                v = np.zeros(m); v[idx[K | Aset]] += 1; v[idx[K | Bset]] += 1; v[idx[K]] -= 1
                gens.append((K, v))
    if EXTRA in ('first', 'sharp'):
        # coarse partitions of K^c into p >= 3 blocks: first-order bound sum_j c_{K+B_j} <= (p-1) c_K + s_K,
        # with Lemma 5.1 s_K <= c_K + sqrt(e e'):  sum_j c_{K+B_j} - p c_K <= sqrt(e_K e'_K).
        # 'sharp' (conjectural, traceless Z): sum_j c_{K+B_j} <= (p-2) c_K + s_K  =>  sum_j c_{K+B_j} - (p-1) c_K <= sqrt.
        def set_partitions(items):
            if not items:
                yield []
                return
            first, rest = items[0], items[1:]
            for part in set_partitions(rest):
                for i in range(len(part)):
                    yield part[:i] + [[first] + part[i]] + part[i + 1:]
                yield [[first]] + part
        for K in subs:
            rest = sorted(full - K)
            if len(rest) < 3:
                continue
            for part in set_partitions(rest):
                if len(part) < 3:
                    continue
                v = np.zeros(m)
                for blk in part:
                    v[idx[K | frozenset(blk)]] += 1
                v[idx[K]] -= len(part) if EXTRA == 'first' else len(part) - 1
                gens.append((K, v))
    G = len(gens)
    # (ii) works modulo P_S = P_{S^c}: map each coordinate S to the class {S, S^c}
    classes = {}
    for S in subs:
        key = min(S, full - S, key=lambda T: (len(T), sorted(T)))
        classes.setdefault(key, []).append(S)
    cl_keys = list(classes); cidx = {}
    for j, key in enumerate(cl_keys):
        for S in classes[key]:
            cidx[S] = j
    mc = len(cl_keys)
    def fold(v):
        out = np.zeros(mc)
        for S in subs:
            out[cidx[S]] += v[idx[S]]
        return out
    # variables: w (G), l (m, for excl_uv), kappa, W-consistency is automatic; m_S (m, excl_ww), kappa'
    nv = G + m + 1 + m + 1
    Aeq, beq = [], []
    # (i): sum_g w_g gen_g - B - sum_R l_R EXCL_R - kappa e_0 = 0  (coordinate-wise)
    for S in subs:
        if S == full:
            continue          # c_[n] = |<u|v>|^2 = 0 for orthogonal u, v: not a free coordinate
        row = np.zeros(nv)
        for g, (K, v) in enumerate(gens):
            row[g] = v[idx[S]]
        for j, R in enumerate(subs):
            row[G + j] = -EXCL[R][idx[S]]
        if len(S) == 0:
            row[G + m] = -1
        Aeq.append(row); beq.append(B[idx[S]])
    # (ii): A(P) - sum_K W_K (P_0 - P_K) - sum_S m_S EXCL_S(P) - kappa' P_0 = 0, folded
    for j in range(mc):
        row = np.zeros(nv)
        for g, (K, v) in enumerate(gens):
            e = np.zeros(m); e[idx[frozenset()]] += 1; e[idx[K]] -= 1     # e_K = P_0 - P_K
            row[g] = fold(e)[j]
        for jj, S in enumerate(subs):
            row[G + m + 1 + jj] = fold(EXCL[S])[j]
        if cl_keys[j] == frozenset():
            row[-1] = 1
        Aeq.append(row); beq.append(fold(A)[j])
    res = linprog(np.zeros(nv), A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=[(0, None)] * nv, method='highs')
    return res.status, res, gens, G
print('extra generators:', EXTRA)
for n in (3, 5):
    status, res, gens, G = run(n)
    print(f'n = {n}: {len(gens)} WZ generators; LP status {status} ({"feasible: certificate exists" if status == 0 else "infeasible"})', flush=True)
    if status == 0:
        w = res.x[:G]
        used = [(sorted(K), round(x, 4)) for (K, v), x in zip(gens, w) if x > 1e-9]
        print('  nonzero WZ multipliers (traced set K, weight):', used[:20], '...' if len(used) > 20 else '')
