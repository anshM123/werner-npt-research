"""Reproduces the failures of sector-restricted (relative) anchors outside the single-excitation sector
(OBSTRUCTIONS.md).  A local U(1) assigns integer weights m(letter); a sector is a span of strings with fixed total weight.
For p in a row sector K, Y = compression of (M - R_p)^Gamma to K (x) K' (rows K, columns K').  We exhibit
(a) same-sector failures in qubit corners: rows and columns both in a full two-letter subspace;
(b) cross-sector failures: p in H_1 with rows in H_2, and a Dicke vector p in H_2 with rows in H_1;
and (c) check that the single-excitation sector itself passes at the same sizes.  Prints PASS if (a), (b) find
negative eigenvalues below -1e-3 and (c) finds none below -1e-12."""
import numpy as np, itertools

def sector(n, j, m): return [s for s in itertools.product(range(3), repeat=n) if sum(m[x] for x in s) == j]

def Mmat(R, C, R2, C2):
    R, C, R2, C2 = map(np.array, (R, C, R2, C2))
    a = R[:, None, None, None, :]; y = C[None, :, None, None, :]; b = R2[None, None, :, None, :]; yp = C2[None, None, None, :, :]
    return np.prod((a == b) * (y == yp) * 1.0 - 0.5 * ((a == y) * (b == yp)), axis=-1)

def lmin(prow, p, rows, cols):
    """p given on the strings prow (its own sector); rows, cols: the compression sets."""
    allc = list(itertools.product(range(3), repeat=len(prow[0])))
    D = np.einsum('a,axbz,b->xz', p.conj(), Mmat(prow, allc, prow, allc), p)
    G = np.einsum('aybx,b->ayx', Mmat(rows, cols, prow, allc), p).reshape(len(rows) * len(cols), len(allc))
    Mc = Mmat(rows, cols, rows, cols).reshape(len(rows) * len(cols), -1)
    W = (Mc - G @ np.linalg.solve(D, G.conj().T)).reshape(len(rows), len(cols), len(rows), len(cols))
    Y = np.transpose(W, (2, 1, 0, 3)).reshape(len(rows) * len(cols), -1)
    return np.linalg.eigvalsh((Y + Y.conj().T) / 2)[0]

rng = np.random.default_rng(0); cz = lambda k: rng.normal(size=k) + 1j * rng.normal(size=k)
diag, two = (0, 1, 1), (0, 0, 1)
res = {}
K = sector(4, 4, diag); res['(a) diagonal charge, n=4, rows H_4, cols H_4'] = min(lmin(K, cz(len(K)), K, K) for _ in range(3))
for n in (3, 4):
    K = sector(n, 0, two); res[f'(a) letter-2 charge, n={n}, rows and cols {{0,1}}^n'] = min(lmin(K, cz(len(K)), K, K) for _ in range(3))
H1, H2 = sector(3, 1, diag), sector(3, 2, diag)
res['(b) n=3, p in H_1 (W), rows H_2, cols H_1'] = lmin(H1, np.array([1.0 if s.count(1) == 1 else 0.0 for s in H1]), H2, H1)
res['(b) n=3, p in H_1 (random), rows H_2, cols H_1'] = min(lmin(H1, cz(len(H1)), H2, H1) for _ in range(3))
res['(b) n=3, p in H_2 (Dicke), rows H_1, cols H_2'] = lmin(H2, np.array([1.0 if s.count(1) == 2 else 0.0 for s in H2]), H1, H2)
ok = True
for k, v in res.items():
    print(f"{k}: smallest eigenvalue {v:+.3e}"); ok &= v < -1e-3
worst = 9
for n in (3, 4):
    H1 = sector(n, 1, diag)
    for c in range(0, n + 1):
        Hc = sector(n, c, diag)
        for _ in range(2): worst = min(worst, lmin(H1, cz(len(H1)), H1, Hc))
print(f"(c) single-excitation sector, n = 3, 4, every column sector, random complex p: smallest {worst:+.2e}"); ok &= worst > -1e-12
print("PASS" if ok else "FAIL")
