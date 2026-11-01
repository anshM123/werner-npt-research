"""Validation of the four-replica block construction.
 (1) true points: for random projectors P, compute Z = twirl(P^(x)4) from the multi-replica moments Tr[P^(x)4 t]
     (t in S_4^n), site-symmetrise, project onto the parameters and check: reconstruction residual ~ 0, every
     constraint block PSD, moments exact, objective = RANK^2 q_n(P) with q_n computed directly from P.
 (2) dense check (n = 1, 2): the twirl is computed independently as the Hilbert-Schmidt projection of P^(x)4 onto the
     span of the 24^n permutation tuples (dense operators on (C^3)^(x)4n); for every constraint the full dense spectrum
     of the constraint operator is compared with the spectrum assembled from the reduced blocks (with multiplicities).
usage: python validate.py n rank [dense] [cons]
"""
import sys, time, itertools, numpy as np
from alg4 import *
from rep4 import Relax, parse_spec, unsvec, tr_embed, orbit_of, constraint_terms

def rand_proj(n, rank, rng, real=False):
    D = 3 ** n
    X = rng.standard_normal((D, rank)) + (0 if real else 1j * rng.standard_normal((D, rank)))
    Q, _ = np.linalg.qr(X)
    return Q @ Q.conj().T

def qn_direct(P, n):
    Pt = P.reshape((3,) * (2 * n)); tot = 0.0
    for k in range(n + 1):
        for J in itertools.combinations(range(n), k):
            rows = [chr(97 + s) for s in range(n)]; cols = [chr(97 + n + s) for s in range(n)]
            for s in J: cols[s] = rows[s]
            out = ''.join(rows[s] for s in range(n) if s not in J) + ''.join(cols[s] for s in range(n) if s not in J)
            R = np.einsum(''.join(rows) + ''.join(cols) + '->' + out, Pt)
            tot += (-0.5) ** k * np.sum(np.abs(R) ** 2)
    return float(tot)

def moments(P, n):
    """m[t] = Tr[P^(x)4 (sigma_1 (x) ... (x) sigma_n)], t = (sigma_s) as S4 indices; site-major convention.
    P is factored as P = U U^dagger (eigenvectors with eigenvalue > 1e-12, scaled) for a cheap tensor-network contraction."""
    ev, V = np.linalg.eigh(P)
    keep = ev > 1e-12
    U = V[:, keep] * np.sqrt(ev[keep])
    k = U.shape[1]
    Ut = U.reshape((3,) * n + (k,)); Uc = Ut.conj()
    out = np.zeros((24,) * n, complex)
    PT_flat = np.array([PERM[s].T.ravel() for s in S4])
    letters = [chr(c) for c in list(range(65, 91)) + list(range(97, 123))]
    a = letters[:4]; b = letters[4:8]; R = letters[8:12]
    X = [[letters[12 + r * n + s] for s in range(n)] for r in range(4)]
    for head in itertools.product(range(24), repeat=n - 1):
        sig = [S4[kk] for kk in head]
        subs = []; ops = []
        for r in range(4):
            subs.append(''.join(X[r][s] for s in range(n - 1)) + a[r] + R[r]); ops.append(Ut)
            subs.append(''.join(X[inv(sig[s])[r]][s] for s in range(n - 1)) + b[r] + R[r]); ops.append(Uc)
        M = np.einsum(','.join(subs) + '->' + ''.join(a) + ''.join(b), *ops, optimize='greedy').reshape(81, 81)
        out[head] = PT_flat @ M.ravel()
    return out

def F1_matrix(A):
    rows = []
    rho = {s: A.coeffs(PERM[s]) for s in S4}
    for lam in range(A.nb):
        f = A.f[lam]
        for i in range(f):
            for j in range(f):
                rows.append([(f / 24.0) * rho[s][lam][j, i] for s in S4])
    return np.array(rows)          # 23 x 24

def zfull_from_moments(m, n, A):
    F1 = F1_matrix(A)
    Z = m
    for s in range(n):
        Z = np.tensordot(F1, Z, axes=([1], [s]))          # new axis 0
        Z = np.moveaxis(Z, 0, s)
    Z = Z.real
    acc = np.zeros_like(Z)
    perms = list(itertools.permutations(range(n)))
    for p in perms: acc += Z.transpose(p)
    return acc / len(perms)

# ---------------------------------------------------------------- dense (n <= 2)
def perm_index(s):
    return np.argmax(PERM[s], axis=0)        # x -> y with PERM[y, x] = 1

def dense_twirl(m, n):
    g1 = np.array([[3.0 ** cycles(mul(inv(a), b)) for b in S4] for a in S4])
    G = np.array([[1.0]])
    for _ in range(n): G = np.kron(G, g1)
    tuples = list(itertools.product(range(24), repeat=n))
    invidx = [S4.index(inv(s)) for s in S4]
    v = np.array([m[tuple(invidx[k] for k in t)] for t in tuples])
    c = np.linalg.pinv(G, rcond=1e-12) @ v
    c = c.real
    D = 81 ** n
    Z = np.zeros((D, D))
    pidx = [perm_index(s) for s in S4]
    for t, ct in zip(tuples, c):
        rows = np.array([0]);
        for k in t: rows = (rows[:, None] * 81 + pidx[k][None, :]).ravel()
        Z[rows, np.arange(D)] += ct
    if n == 2:
        Z = (Z + Z.reshape(81, 81, 81, 81).transpose(1, 0, 3, 2).reshape(D, D)) / 2
    return Z

def dense_map(Z, n, T, ST):
    """at every site: trace the replicas in T (identities put back), then partial transpose on the replicas in ST"""
    D = 81 ** n
    X = Z.reshape((3,) * (8 * n))
    letters = [chr(c) for c in list(range(97, 123)) + list(range(65, 91))]
    row = letters[:4 * n]; col = letters[4 * n:8 * n]
    if T:
        colt = list(col)
        for s in range(n):
            for r in T: colt[4 * s + r] = row[4 * s + r]
        keep = [k for k in range(4 * n) if (k % 4) not in T]
        out = ''.join(row[k] for k in keep) + ''.join(col[k] for k in keep)
        Y = np.einsum(''.join(row) + ''.join(colt) + '->' + out, X)
        subs = [out]; ops = [Y]
        for s in range(n):
            for r in T: subs.append(row[4 * s + r] + col[4 * s + r]); ops.append(np.eye(3))
        X = np.einsum(','.join(subs) + '->' + ''.join(row) + ''.join(col), *ops)
    if ST:
        perm = list(range(8 * n))
        for s in range(n):
            for r in ST: perm[4 * s + r], perm[4 * n + 4 * s + r] = perm[4 * n + 4 * s + r], perm[4 * s + r]
        X = X.transpose(perm)
    return np.ascontiguousarray(X).reshape(D, D)


def dense_leg_perm(Z, n, p):
    """apply the per-site leg permutation p (new leg k = old leg p[k]; legs 0-3 rows, 4-7 cols) at every site"""
    D = 81 ** n
    X = Z.reshape((3,) * (8 * n))
    def ax(s, k): return 4 * s + k if k < 4 else 4 * n + 4 * s + (k - 4)
    perm = [0] * (8 * n)
    for s in range(n):
        for k in range(8): perm[ax(s, k)] = ax(s, p[k])
    return np.ascontiguousarray(X.transpose(perm)).reshape(D, D)

def dense_constraint(Zd, n, cname, RANK):
    SQ, ST, R = parse_spec(cname)
    tgkey, tdefs = constraint_terms(cname, RANK)
    barred = set(tgkey)
    if R == 0: b0 = set(ST)
    elif R == 1: b0 = set(ST) | {3}
    else: b0 = {2, 3}
    order = [r for r in range(4) if r not in b0] + sorted(b0)
    L = np.zeros_like(Zd)
    for k in range(len(SQ) + 1):
        for T in itertools.combinations(SQ, k):
            Y = dense_map(Zd, n, list(T), list(ST))
            if R == 1: Y = dense_leg_perm(Y, n, (0, 1, 2, 6, 4, 5, 7, 3))
            elif R == 2: Y = dense_leg_perm(Y, n, (0, 2, 4, 6, 5, 7, 1, 3))
            if order != [0, 1, 2, 3]: Y = dense_leg_perm(Y, n, tuple(order) + tuple(4 + o for o in order))
            L += (-1.0) ** (len(SQ) - k) / float(RANK) ** k * Y
    return L, tgkey

def block_spectrum(rel, w, cname, tg):
    ev = []
    if cname == 'Z':
        for s in rel.src:
            mult = np.prod([rel.A.m[l] for l in s['lam']]) * len(s['orbit'])
            k0 = 0
            for c in s['classes']:
                k = c['m'] * (c['m'] + 1) // 2
                W = unsvec(w[s['off'] + k0: s['off'] + k0 + k], c['m']); k0 += k
                for e in np.linalg.eigvalsh(W): ev += [e / np.prod([rel.A.m[l] for l in s['lam']])] * int(mult * c['d'])
        return np.sort(np.array(ev))
    for (cn, mu, c, M) in rel.blocks:
        if cn != cname: continue
        mt = np.prod([tg.m[l] for l in mu]); mult = mt * len(orbit_of(mu)) * c['d']
        V = unsvec(M @ w, c['m'])
        for e in np.linalg.eigvalsh(V): ev += [e / mt] * int(mult)
    ev = np.array(ev)
    D = 81 ** rel.n
    ev = np.concatenate([ev, np.zeros(D - len(ev))])
    return np.sort(ev)

if __name__ == "__main__":
    n = int(sys.argv[1]); RANK = int(sys.argv[2])
    dense = len(sys.argv) > 3 and sys.argv[3] == 'dense'
    cons = sys.argv[4].split(',') if len(sys.argv) > 4 else ['LE', 'PT1', 'PT2', 'PT1LE', 'PT2LE', 'QQ', 'QQQ', 'QQQQ', 'PROJ']
    algs = build_algebras(); A = algs[0]
    t0 = time.time()
    rel = Relax(n, RANK, algs, cons, verbose=False)
    print(f"n={n} rank={RANK}: relaxation built ({rel.N} params) in {time.time() - t0:.1f}s; constraints {cons}", flush=True)
    rng = np.random.default_rng(100 + n)
    for trial in range(2 if not dense else 1):
        P = rand_proj(n, RANK, rng, real=(trial == 1))
        t1 = time.time()
        m = moments(P, n)
        Zf = zfull_from_moments(m, n, A)
        w, res = rel.w_from_zfull(Zf)
        mins, mres, eres, objv = rel.evaluate(w)
        q = qn_direct(P, n)
        print(f"  trial {trial} ({'real' if trial == 1 else 'complex'} P): reconstruction residual {res:.1e}; "
              f"moment residual {mres:.1e}; PROJ residual {eres:.1e}; objective/R^2 = {objv / RANK ** 2:.12f}, "
              f"direct q_n(P) = {q:.12f} (diff {abs(objv / RANK ** 2 - q):.1e}); {time.time() - t1:.1f}s", flush=True)
        print("     min eigenvalue per constraint at the true point:", {k: f"{v:.2e}" for k, v in mins.items()}, flush=True)
        if dense:
            t2 = time.time()
            Zd = dense_twirl(m, n)
            print(f"  dense twirl built ({Zd.shape[0]} dim) in {time.time() - t2:.1f}s", flush=True)
            ed = np.sort(np.linalg.eigvalsh(Zd)); eb = block_spectrum(rel, w, 'Z', A)
            print(f"     Z: dense vs block spectrum max diff {np.abs(ed - eb).max():.1e}  (dense min {ed[0]:.3e})", flush=True)
            if 'PROJ' in cons:
                F = PERM[(0, 1, 3, 2)]
                if n == 1:
                    Lp = ptrace_last(Zd @ F) - ptrace_last(Zd) / RANK
                else:
                    Ff = np.kron(F, F)
                    X = (Zd @ Ff).reshape(27, 3, 27, 3, 27, 3, 27, 3)
                    Y = Zd.reshape(27, 3, 27, 3, 27, 3, 27, 3)
                    Lp = np.einsum('aibjcidj->abcd', X) - np.einsum('aibjcidj->abcd', Y) / RANK
                print(f"     PROJ: dense identity residual {np.abs(Lp).max():.1e}", flush=True)
            for cname in cons:
                if cname == 'PROJ': continue
                L, tgkey = dense_constraint(Zd, n, cname, RANK)
                tg = {(): algs[0], (3,): algs[1], (2, 3): algs[2]}[tgkey]
                ed = np.sort(np.linalg.eigvalsh((L + L.T) / 2)); eb = block_spectrum(rel, w, cname, tg)
                print(f"     {cname}: dense vs block spectrum max diff {np.abs(ed - eb).max():.1e} (dense min eig {ed[0]:.3e}, "
                      f"asym {np.abs(L - L.T).max():.1e}); {time.time() - t2:.0f}s", flush=True)
