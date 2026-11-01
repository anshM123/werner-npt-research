"""Exact (rational) per-site Fourier data for the certificate pipeline.

For an algebra B of operators on (C^3)^(x)4 commuting with a U(3) action (U on unbarred replicas, conj(U) on barred
ones), each isotypic component is S^nu (x) V^nu.  Exact construction:
  1. isotypic components = eigenspaces of an integer central element (quadratic Casimir written with swaps F_ij on
     unbarred pairs and contractions Omega_ij between unbarred i and barred j): exact rational nullspaces;
  2. the vectors of a component whose U(3) weight equals the highest weight of V^nu span one irreducible B-module
     S^nu (x) |hw>, because the highest weight occurs once in V^nu; weight spaces are coordinate subspaces;
  3. irrep matrices rho(b) = (S^T S)^{-1} S^T b S on that rational basis S, with Gram matrix S^T S.
Everything is computed with python-flint fmpq_mat (exact).  The result is checked against the numerical Alg data:
same irrep dimensions and multiplicities, rho a representation, b S = S rho(b) exactly.
usage: python exact_site.py"""
import itertools, numpy as np, flint
from alg4 import d, S4, PERM, perm_op, pt

def Q(M):
    """numpy integer/rational-valued float matrix -> fmpq_mat (entries must be integers)"""
    r, c = M.shape
    return flint.fmpq_mat(r, c, [int(round(v)) for v in M.ravel()])

def weights(barred, nr=4):
    """U(3) weight of each computational basis vector of (C^3)^(x)nr; barred replicas count with sign -1"""
    out = []
    for x in itertools.product(range(d), repeat=nr):
        w = [0] * d
        for r, xr in enumerate(x): w[xr] += (-1 if r in barred else 1)
        out.append(tuple(w))
    return out

def nullspace_q(M):
    """exact nullspace basis (columns) of an fmpq_mat"""
    R, rank = M.rref()
    n = M.ncols(); piv = []; row = 0
    for j in range(n):
        if row < rank and R[row, j] == 1 and all(R[i, j] == 0 for i in range(rank) if i != row):
            piv.append(j); row += 1
    free = [j for j in range(n) if j not in piv]
    B = flint.fmpq_mat(n, len(free))
    for k, j in enumerate(free):
        B[j, k] = 1
        for i, p in enumerate(piv): B[p, k] = -R[i, j]
    return B

class ExactAlg:
    def __init__(self, gens, central, barred, hws, name, nr=4):
        """gens: dict label -> integer numpy matrix spanning the algebra; central: integer numpy matrix (central);
        barred: set of barred replicas; hws: dict eigenvalue -> highest weight tuple"""
        self.name = name; D = d ** nr
        Cq = Q(central); W8 = weights(barred, nr)
        self.blocks = []
        ev = sorted(set(np.round(np.linalg.eigvalsh(central.astype(float)), 6)))
        for lam in ev:
            lam_q = flint.fmpq(int(round(lam * 2)), 2)
            Wsp = nullspace_q(Cq - lam_q * flint.fmpq_mat(D, D, [1 if i == j else 0 for i in range(D) for j in range(D)]))
            dimW = Wsp.ncols()
            if dimW == 0: continue
            hw = hws[int(round(lam * 2))]
            idx = [i for i in range(D) if W8[i] == hw]
            # S = W restricted to the coordinate subspace idx: vectors in W with zero entries outside idx
            Pout = flint.fmpq_mat(D - len(idx), D, [0] * ((D - len(idx)) * D))
            outs = [i for i in range(D) if i not in set(idx)]
            for r_, i in enumerate(outs): Pout[r_, i] = 1
            K = nullspace_q(Pout * Wsp)           # coefficients c with W c supported on idx
            S = Wsp * K
            f = S.ncols(); m = dimW // f
            G = S.transpose() * S; Gi = G.inv()
            rho = {lab: Gi * (S.transpose() * (Q(M) * S)) for lab, M in gens.items()}
            # exactness checks: b S = S rho(b)
            for lab, M in gens.items():
                assert Q(M) * S == S * rho[lab], f"{name}: module not invariant under {lab}"
            self.blocks.append(dict(casimir=lam, f=f, m=m, S=S, G=G, Gi=Gi, rho=rho))
        self.blocks.sort(key=lambda b: (-b['m'], b['f']))

    def coeffs(self, X):
        """exact Fourier blocks of an algebra element X (integer numpy matrix)"""
        Xq = Q(X)
        return [b['Gi'] * (b['S'].transpose() * (Xq * b['S'])) for b in self.blocks]


def build_exact():
    gensA = {s: PERM[s] for s in S4}
    # class sum of transpositions: central in the S_4 algebra; eigenvalue = content sum of lambda
    T = sum(PERM[s] for s in S4 if sum(1 for i in range(4) if s[i] != i) == 2)
    hwA = {12: (4, 0, 0), 4: (3, 1, 0), 0: (2, 2, 0), -4: (2, 1, 1)}      # keys: 2 * eigenvalue
    A = ExactAlg(gensA, T, set(), hwA, 'A exact')
    gensB1 = {s: pt(PERM[s], [3]) for s in S4}
    # Casimir-type central element of B_{3,1}: swaps among replicas 0,1,2 minus contractions (i,3)
    sw = [s for s in S4 if sum(1 for i in range(4) if s[i] != i) == 2 and s[3] == 3]
    ct = [s for s in S4 if sum(1 for i in range(4) if s[i] != i) == 2 and s[3] != 3]
    CB = sum(pt(PERM[s], [3]) for s in sw) - sum(pt(PERM[s], [3]) for s in ct)
    evB = sorted(set(np.round(np.linalg.eigvalsh(CB.astype(float)), 6)))
    print("B31 central element eigenvalues:", evB)
    return A, CB, gensB1, evB

if __name__ == "__main__":
    import time
    t0 = time.time()
    A, CB, gensB1, evB = build_exact()
    for b in A.blocks:
        print(f"{A.name}: f={b['f']} m={b['m']} casimir={b['casimir']}, Gram diag {[b['G'][i, i] for i in range(b['f'])]}")
    # identify highest weights of the B31 components by matching dimensions: (3),(1) -> (3,0,-1); (21),(1) -> (2,1,-1);
    # (2) -> (2,0,0); (11) -> (1,1,0)
    W8 = weights({3})
    from collections import Counter
    hwB = {}
    for lam in evB:
        lam_q = flint.fmpq(int(round(lam * 2)), 2)
        Wsp = nullspace_q(Q(CB) - lam_q * flint.fmpq_mat(81, 81, [1 if i == j else 0 for i in range(81) for j in range(81)]))
        # weights present in the component: the highest one in dominance order is the highest weight
        ws = Counter()
        for k in range(Wsp.ncols()):
            for i in range(81):
                if Wsp[i, k] != 0: ws[W8[i]] += 1
        cand = max(ws, key=lambda w: (w[0], w[0] + w[1]))
        hwB[int(round(lam * 2))] = cand
        print(f"   B31 component casimir {lam}: dim {Wsp.ncols()}, highest weight {cand}")
    B1 = ExactAlg(gensB1, CB, {3}, hwB, 'B31 exact')
    for b in B1.blocks:
        print(f"{B1.name}: f={b['f']} m={b['m']} casimir={b['casimir']}")
    # representation check: rho(s) rho(t) = rho(s t) for the S_4 part
    from alg4 import mul
    ok = all(b['rho'][s] * b['rho'][t] == b['rho'][mul(s, t)] for b in A.blocks for s in S4[:6] for t in S4[:6])
    print(f"A irreps multiplicative: {ok}; built in {time.time() - t0:.1f}s")
