"""A check of the marked-copy theorems (letter-1 block of every (R_inf)_j at weights (a, 1, ..., 1), sigma = a^2),
written separately from the certificates: it uses only the Lemma B builder j1lib.py (for its operators and, at small
k, for an entry-by-entry comparison), python-flint and numpy.  None of the reduction (twisted invariants, orbit
calculus, families) is used.

  python markcheck.py V K A1,A2,...     the fast builder below against J1.matrix (Lemma B, entry by entry), k = K
  python markcheck.py X K A1,A2,...     exact: for every j = -1..K and every a, the known kernel vectors are
                                        annihilated (exact), the block without a set of pivot coordinates of them is
                                        positive definite (pd_cert, exact), and the number of known kernel vectors is
                                        C(K, j) + [j = K-2] + (K-1)[j = K-1]
  python markcheck.py N K A1,A2,...     negative control: B scaled by 1001/1000 (and, for K <= 5, the scaled builder
                                        against J1 with bscale); at least one complement must be refused by pd_cert
  python markcheck.py F K A1,A2,...     floating point (entries exact, then float64): kernel vectors annihilated
                                        exactly; eigenvalues of the block: the number below 1e-9 * max equals the
                                        number of known kernel vectors and none is below -1e-9 * max (not a proof)

Fast builder: the B-part of Lemma B is a Gram form of E^{-1} on the vectors T_r|y> and J^+|y> (E symmetric,
J^- = (J^+)^T), computed sector by sector with flint matrices; the index pattern is that of J1.entry.
"""
import sys, time, itertools, math
from fractions import Fraction as Fr
from math import comb
import numpy as np
import flint

from crosscheck import pd_cert, pivots_of, known_kernel
from j1lib import J1


def log(*a):
    print(*a, flush=True)


def ones_strings(k, c):
    out = []
    for P in itertools.combinations(range(k), c):
        y = [0] * k
        for t in P:
            y[t] = 1
        out.append(tuple(y))
    return out


def fq(x):
    x = Fr(x)
    return flint.fmpq(x.numerator, x.denominator)


def tofr(x):
    return Fr(int(x.p), int(x.q))


class Fast:
    def __init__(self, w, bscale=1):
        self.bscale = Fr(bscale)
        self.w = [Fr(x) for x in w]
        self.k = len(w)
        self.J = J1(self.w)
        self.S = {}
        self.idx = {}
        self.Einv = {}

    def strings(self, c):
        if c not in self.S:
            self.S[c] = ones_strings(self.k, c) if 0 <= c <= self.k else []
            self.idx[c] = {y: i for i, y in enumerate(self.S[c])}
        return self.S[c]

    def E(self, c):
        """E on the strings with c letters 1, as a flint fmpq_mat (from J1.Eop)"""
        S = self.strings(c)
        idx = self.idx[c]
        n = len(S)
        ent = [[flint.fmpq(0)] * n for _ in range(n)]
        for jj, y in enumerate(S):
            for z, cf in self.J.Eop({y: Fr(1)}).items():
                ent[idx[z]][jj] += fq(cf)
        return flint.fmpq_mat(n, n, [x for row in ent for x in row])

    def cols(self, c, vecs):
        """matrix with columns the dict vectors (all supported on strings with c letters 1)"""
        S = self.strings(c)
        idx = self.idx[c]
        n, m = len(S), len(vecs)
        ent = [[flint.fmpq(0)] * m for _ in range(n)]
        for jj, v in enumerate(vecs):
            for z, cf in v.items():
                assert sum(1 for t in z if t == 1) == c and all(t in (0, 1) for t in z)
                ent[idx[z]][jj] += fq(cf)
        return flint.fmpq_mat(n, m, [x for row in ent for x in row])

    def gram(self, c, P, Q):
        """P^T E_c^{-1} Q for column matrices P, Q on the c-strings"""
        if P.ncols() == 0 or Q.ncols() == 0 or P.nrows() == 0:
            return flint.fmpq_mat(P.ncols(), Q.ncols())
        E = self.E(c)
        return P.transpose() * E.solve(Q)

    def block(self, j):
        """letter-1 group of (R_inf)_j: states ((1, r), y), |y| = j, then (None, y), |y| = j + 1; exact Fractions"""
        k, J = self.k, self.J
        S1 = self.strings(j) if 0 <= j <= k else []
        S2 = self.strings(j + 1) if 0 <= j + 1 <= k else []
        self.strings(j + 2)
        st1 = [((1, r), y) for y in S1 for r in range(k)]
        st0 = [(None, y) for y in S2]
        states = st1 + st0
        n1, n0 = len(st1), len(st0)
        # G11[(r1, y1), (r2, y2)] = <T_r1 y1, E_j^{-1} T_r2 y2>
        if n1:
            X1 = self.cols(j, [J.Tm((1, r), {y: Fr(1)}) for (_, r), y in st1])
            G11 = self.gram(j, X1, X1)
        # G12[y', (r, y)] = <J^+ y', E_{j+1}^{-1} T_r y>,  y' in S1, y in S2
        if n1 and n0:
            pairs2 = [(r, y) for y in S2 for r in range(k)]
            X2 = self.cols(j + 1, [J.Tm((1, r), {y: Fr(1)}) for r, y in pairs2])
            P1 = self.cols(j + 1, [J.Jp({y: Fr(1)}) for y in S1])
            G12 = self.gram(j + 1, P1, X2)
            pidx2 = {p: i for i, p in enumerate(pairs2)}
        # G22[y', y] = <J^+ y', E_{j+2}^{-1} J^+ y>
        if n0:
            if j + 2 <= k:
                P2 = self.cols(j + 2, [J.Jp({y: Fr(1)}) for y in S2])
                G22 = self.gram(j + 2, P2, P2)
            else:
                G22 = None
        i1 = {s: i for i, s in enumerate(st1)}
        s1idx = {y: i for i, y in enumerate(S1)}
        M = [[None] * (n1 + n0) for _ in range(n1 + n0)]
        for p, (a, yp) in enumerate(st1):
            ra = a[1]
            for q, (b, y) in enumerate(st1):
                if q < p:
                    continue
                rb = b[1]
                coef = Fr(1, 2) - (1 if ra == rb else 0)
                loc = Fr(0)
                if yp[rb] == 1 and y[ra] == 1:
                    z1 = list(yp); z1[rb] = 0
                    z2 = list(y); z2[ra] = 0
                    if z1 == z2:
                        loc += coef
                if a == b and yp == y:
                    loc += 2 - (1 if y[ra] != 0 else 0)
                Bv = self.bscale * tofr(G11[i1[(b, yp)], i1[(a, y)]])
                M[p][q] = M[q][p] = loc - Bv
            for q0, (_, y) in enumerate(st0):
                q = n1 + q0
                sg = 0
                if y[ra] == 1:
                    z = list(y); z[ra] = 0
                    if tuple(z) == yp:
                        sg = 1
                Bv = self.bscale * tofr(G12[s1idx[yp], pidx2[(ra, y)]])
                M[p][q] = M[q][p] = -sg + Bv
        for p0, (_, yp) in enumerate(st0):
            for q0, (_, y) in enumerate(st0):
                if q0 < p0:
                    continue
                Bv = self.bscale * tofr(G22[p0, q0]) if G22 is not None else Fr(0)
                M[n1 + p0][n1 + q0] = M[n1 + q0][n1 + p0] = 2 * (1 if yp == y else 0) - Bv
        return states, M


def kernel_ok(M, vecs):
    n = len(M)
    if not vecs:
        return True
    A = flint.fmpq_mat(n, n, [fq(x) for row in M for x in row])
    V = flint.fmpq_mat(n, len(vecs), [fq(vecs[c][r]) for r in range(n) for c in range(len(vecs))])
    P = A * V
    return all(P[i, c] == 0 for i in range(n) for c in range(len(vecs)))


def expected_nker(k, j):
    return (comb(k, j) if 0 <= j <= k else 0) + (1 if j == k - 2 else 0) + (k - 1 if j == k - 1 else 0)


def main():
    mode, k = sys.argv[1], int(sys.argv[2])
    avals = [Fr(x) for x in sys.argv[3].split(',')]
    allok = True
    nrefused = 0
    for a in avals:
        w = [a] + [Fr(1)] * (k - 1)
        fb = Fast(w, Fr(1001, 1000) if mode == 'N' else 1)
        for j in range(-1, k + 1):
            t0 = time.time()
            states, M = fb.block(j)
            n = len(states)
            if mode == 'V':
                Ms = J1(w).matrix(states)
                ok = all(M[p][q] == Ms[p][q] for p in range(n) for q in range(n))
                log(f'V k={k} a={a} j={j}: {n} states, fast builder = J1.matrix entry by entry: {ok} [{time.time() - t0:.1f}s]')
                allok &= ok
                continue
            vecs = known_kernel(states, k, j, w)
            kern = kernel_ok(M, vecs)
            cnt = len(vecs) == expected_nker(k, j)
            piv = set(pivots_of(vecs)) if vecs else set()
            keep = [i for i in range(n) if i not in piv]
            if mode == 'N':
                Ms = J1(w, bscale=Fr(1001, 1000)).matrix(states) if k <= 5 else None
                same = Ms is None or all(M[p][q] == Ms[p][q] for p in range(n) for q in range(n))
                proved, info = pd_cert([[M[i][jj] for jj in keep] for i in keep])
                nrefused += (not proved)
                allok &= same
                log(f'N k={k} a={a} (sigma={a * a}) j={j}: B scaled by 1001/1000 (builder = J1 with bscale: {same}); '
                    f'complement certified PD {proved} ({info})')
                continue
            if mode == 'X':
                proved, info = pd_cert([[M[i][jj] for jj in keep] for i in keep])
                ok = kern and cnt and proved
                log(f'X k={k} a={a} (sigma={a * a}) j={j}: {n} states, {len(vecs)} known kernel vectors (expected '
                    f'{expected_nker(k, j)}), annihilated {kern}, PD without pivots {proved} ({info}) '
                    f'[{time.time() - t0:.1f}s]')
            else:
                Af = np.array([[float(x) for x in row] for row in M]) if n else np.zeros((0, 0))
                ev = np.linalg.eigvalsh(Af) if n else np.zeros(0)
                mx = max(1.0, float(np.max(np.abs(ev)))) if n else 1.0
                nz = int(np.sum(np.abs(ev) < 1e-9 * mx))
                neg = int(np.sum(ev < -1e-9 * mx))
                pos = ev[ev > 1e-9 * mx]
                gap = float(pos[0]) if len(pos) else float('nan')
                ok = kern and cnt and nz == len(vecs) and neg == 0
                log(f'F k={k} a={a} (sigma={a * a}) j={j}: {n} states, {len(vecs)} known kernel vectors, annihilated '
                    f'{kern}; eigenvalues: {nz} near 0, {neg} negative, smallest positive {gap:.3e}, largest {mx:.3e}: '
                    f'{ok} [{time.time() - t0:.1f}s]')
            allok &= ok
    if mode == 'N':
        log(f'N: scaled blocks refused for {nrefused} (j, a) pairs (at least one needed)')
        allok &= nrefused >= 1
    log(f'total {time.time() - T0:.0f}s')
    log('PASS' if allok else 'FAIL')


T0 = time.time()
if __name__ == '__main__':
    main()
