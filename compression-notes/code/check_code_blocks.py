"""Code blocks for the qutrit Werner endpoint (rho = I - F/2 on each pair of copies).

For a qutrit stabiliser group S on n sites (|S| = 3^(n-k), code space Q of dimension 3^k), let
sigma_Q = (P_Q (x) I) rho^{(x)n} (P_Q (x) I) be the one-sided compression.
  (a) Weight formula: sigma_Q is block diagonal over the syndrome sectors of the second side, and in each sector its
      partial transpose is diagonal in a basis of maximally entangled vectors labelled by the error classes E + S
      of that syndrome, with eigenvalue p(E + S) = |S|^{-1} sum_{s in S} (-1/2)^{n - wt(E + s)}.
      Checked here against the dense partial transpose for [[3,1,2]] and [[4,2,2]].
  (b) Separability: uniform weights on the 3^k classes of E0 + L (L a maximal isotropic subspace containing S) are
      the twirls of stabiliser product states, hence separable. An LP finds a nonnegative combination of these
      equal to p in every sector; the solution is then made exact (sympy) and re-checked with Fractions.
  (c) The family [[n, n-2, 2]] (stabilisers X^n and Z^c, all c_i nonzero, sum c_i = 0 mod 3): weights depend only on
      the type of E, and the minimum weight is 1 - n/6 for n <= 9, so the compression is NPT for n >= 7.
Paulis mod phase are vectors (a|b) in F_3^(2n)."""
import itertools, numpy as np, sympy as sp
from fractions import Fraction
from scipy.optimize import linprog
from scipy.sparse import csc_matrix


def symp(u, v, n):
    return int((u[:n] @ v[n:] - u[n:] @ v[:n]) % 3)


class Code:
    def __init__(self, gens, n):
        self.n = n
        self.gens = [np.array(g) % 3 for g in gens]
        self.k = n - len(gens)
        S = [np.zeros(2 * n, int)]
        for g in self.gens:
            S = [(s + c * g) % 3 for s in S for c in range(3)]
        self.S = np.array(S)
        assert len({tuple(s) for s in S}) == 3 ** (n - self.k), "generators dependent"
        assert all(symp(g, h, n) == 0 for g in self.gens for h in self.gens), "generators do not commute"
        self.pw = 3 ** np.arange(2 * n)[::-1]

    def enc(self, V):
        return (V % 3) @ self.pw

    def cls(self, v):
        return int(self.enc((v[None, :] + self.S) % 3).min())

    def wt(self, V):
        V = np.atleast_2d(V)
        n = self.n
        return ((V[:, :n] != 0) | (V[:, n:] != 0)).sum(1)

    def weight(self, v):
        w = self.wt((v[None, :] + self.S) % 3)
        return sum((Fraction(-1, 2) ** int(self.n - x) for x in w), Fraction(0)) / len(self.S)

    def syndrome(self, v):
        return tuple(symp(v, g, self.n) for g in self.gens)


def dec(x, n):
    d = []
    for _ in range(2 * n):
        d.append(x % 3)
        x //= 3
    return np.array(d[::-1])


def sectors(code):
    n = code.n
    seen = set()
    out = {}
    for x in range(9 ** n):
        v = dec(x, n)
        c = code.cls(v)
        if c in seen:
            continue
        seen.add(c)
        out.setdefault(code.syndrome(v), []).append(c)
    return out


def logical_lagrangians(code):
    """maximal isotropic subspaces L of F_3^(2n) containing S, as sorted lists of class codes (L / S)."""
    n, k = code.n, code.k
    reps = {}
    for x in range(9 ** n):
        v = dec(x, n)
        if all(symp(v, g, n) == 0 for g in code.gens):
            reps.setdefault(code.cls(v), v)
    reps = list(reps.items())                        # the 9^k logical classes
    zero = code.cls(np.zeros(2 * n, int))
    level = {frozenset([zero]): []}
    for _ in range(k):
        nxt = {}
        for sub, vs in level.items():
            for c, v in reps:
                if c in sub or any(symp(v, w, n) for w in vs):
                    continue
                nv = vs + [v]
                spn = frozenset(code.cls(sum(a * w for a, w in zip(t, nv)) % 3)
                                for t in itertools.product(range(3), repeat=len(nv)))
                if spn not in nxt:
                    nxt[spn] = nv
        level = nxt
    return [sorted(sub) for sub in level]


def certify(code, name, exact=True):
    n, k = code.n, code.k
    secs = sectors(code)
    lags = logical_lagrangians(code)
    print(f"{name}: n = {n}, k = {k}, {len(secs)} sectors of {3 ** (2 * k)} classes, "
          f"{len(lags)} maximal isotropic subspaces containing S")
    allok = True
    wmin = None
    values = set()
    for s, cl in sorted(secs.items()):
        reps = {c: dec(c, n) for c in cl}
        p = {c: code.weight(reps[c]) for c in cl}
        values |= set(p.values())
        wmin = min(p.values()) if wmin is None else min(wmin, min(p.values()))
        idx = {c: i for i, c in enumerate(cl)}
        cols = set()
        for sub in lags:
            for c in cl:
                cols.add(frozenset(code.cls((reps[c] + dec(l, n)) % 3) for l in sub))
        cols = sorted(cols, key=sorted)
        rows, cc = [], []
        for j, col in enumerate(cols):
            for c in col:
                rows.append(idx[c])
                cc.append(j)
        A = csc_matrix((np.full(len(rows), 1.0 / 3 ** k), (rows, cc)), shape=(len(cl), len(cols)))
        b = np.array([float(p[c]) for c in cl])
        res = linprog(np.zeros(len(cols)), A_eq=A, b_eq=b, bounds=(0, None), method="highs-ds")
        if res.status != 0:
            allok = False
            print(f"  sector {s}: p in [{min(p.values())}, {max(p.values())}]; LP infeasible "
                  f"(stabiliser cosets do not suffice)")
            continue
        if not exact:
            print(f"  sector {s}: p in [{min(p.values())}, {max(p.values())}]; LP feasible")
            continue
        supp = [j for j in range(len(cols)) if res.x[j] > 1e-10]
        As = sp.Matrix([[sp.Rational(1, 3 ** k) if c in cols[j] else 0 for j in supp] for c in cl])
        bs = sp.Matrix([sp.Rational(p[c].numerator, p[c].denominator) for c in cl])
        sol, params = As.gauss_jordan_solve(bs)
        sol = sol.subs({t: 0 for t in params})
        lam = [Fraction(int(x.p), int(x.q)) for x in sol]
        tot = {c: Fraction(0) for c in cl}
        for j, l in zip(supp, lam):
            for c in cols[j]:
                tot[c] += l / 3 ** k
        ok = all(l >= 0 for l in lam) and all(tot[c] == p[c] for c in cl)
        allok &= ok
        print(f"  sector {s}: p in [{min(p.values())}, {max(p.values())}]; exact certificate with "
              f"{sum(1 for l in lam if l)} cosets: {'VALID' if ok else 'FAILED'}")
    print(f"  weights {[str(v) for v in sorted(values)]}; minimum {wmin}; "
          f"separable (exact certificates in every sector): {allok}")
    return allok


def dense_check(code, name):
    """dense one-sided compression; PT spectrum against the weight formula (as multisets)."""
    n = code.n
    w = np.exp(2j * np.pi / 3)
    X = np.roll(np.eye(3), 1, axis=0)
    Z = np.diag([1, w, w ** 2])

    def pauli(v):
        M = np.eye(1)
        for i in range(n):
            M = np.kron(M, np.linalg.matrix_power(X, v[i]) @ np.linalg.matrix_power(Z, v[n + i]))
        return M
    D = 3 ** n
    P = np.eye(D)
    for g in code.gens:
        G = pauli(g)
        P = P @ (np.eye(D) + G + G @ G) / 3
    ev, U = np.linalg.eigh((P + P.conj().T) / 2)
    V = U[:, ev > 0.5]
    kk = V.shape[1]
    F = np.zeros((9, 9))
    for i in range(3):
        for j in range(3):
            F[3 * j + i, 3 * i + j] = 1
    rho = np.eye(9) - F / 2
    R = np.eye(1)
    for _ in range(n):
        R = np.kron(R, rho)                                # ordering (A1 B1)(A2 B2)...
    perm = [2 * i for i in range(n)] + [2 * i + 1 for i in range(n)]
    R = R.reshape([3] * (4 * n)).transpose(perm + [2 * n + q for q in perm]).reshape(D * D, D * D)
    E = np.kron(V, np.eye(D))
    sig = E.conj().T @ R @ E
    sg = np.linalg.eigvalsh(sig.reshape(kk, D, kk, D).transpose(0, 3, 2, 1).reshape(kk * D, kk * D))
    form = []
    for s, cl in sectors(code).items():
        for c in cl:
            form.append(float(code.weight(dec(c, n))))
    form = np.sort(form)
    print(f"{name}: dense PT spectrum vs weight formula: max difference {np.abs(np.sort(sg) - form).max():.1e} "
          f"({len(form)} eigenvalues)")


if __name__ == "__main__":
    def Xall(n):
        v = np.zeros(2 * n, int)
        v[:n] = 1
        return v

    def Zc(n, c):
        v = np.zeros(2 * n, int)
        v[n:] = c
        return v
    c312 = Code([Xall(3), Zc(3, [1, 1, 1])], 3)
    c422 = Code([Xall(4), Zc(4, [1, 2, 1, 2])], 4)
    print("(a) weight formula")
    dense_check(c312, "[[3,1,2]]")
    dense_check(c422, "[[4,2,2]]")
    print("(b) separability certificates")
    certify(c312, "[[3,1,2]]")
    certify(c422, "[[4,2,2]]")
    print("(c) [[n, n-2, 2]] family: minimum over types of sum_j (-1/2)^{N_j} / 9")
    for n in range(2, 10):
        types = {tuple(sorted(t)) for t in itertools.combinations_with_replacement(range(n + 1), 9) if sum(t) == n}
        best = min(sum(Fraction(-1, 2) ** Nj for Nj in N) / 9 for N in types)
        print(f"  n = {n}: minimum weight {best}; equals 1 - n/6: {best == 1 - Fraction(n, 6)}")
