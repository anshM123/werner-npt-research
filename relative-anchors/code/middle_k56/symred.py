"""Isotypic reduction of a torus block under the permutations of the copies that preserve the levels and the set T.

G = prod_i Sym(B_i) for disjoint blocks B_i of copies (in configuration (a, b): the copies of T on level 1, of T on
level q, of U on level 1, of U on level q).  An irrep of G is a tuple of partitions lam_i of |B_i|.  Fill the Young
diagram of lam_i row by row with the copies of B_i in increasing order; R and C are the row and column groups
(products over i), and e = (sum_{r in R} r)(sum_{c in C} sgn(c) c) is the Young symmetrizer.  On every permutation
module V of G, e V has dimension mult(lam, V), and its vectors are s0 (x) c for a fixed vector s0 of the irrep
(von Neumann's lemma: e x e is a multiple of e for every x in the group algebra).  Irreps of symmetric groups are
absolutely irreducible and self-dual, so a G-invariant symmetric form Y acts on the lam-isotypic part as <,> (x) Q_lam,
and for a G-stable subspace W and any basis B of e W, B Y B^T = |s0|^2 C Q_lam C^T with C invertible.  Hence
Y|_W > 0 iff B Y B^T > 0 for every lam (the dimension count sum_lam dim(lam) dim(e W) = dim W shows no irrep is missed).

The permutation pi acts on strings by (pi y)[pi(t)] = y[t], on states by acting on row and column, and on vectors by
pi e_x = e_{pi x}."""
import itertools
from fractions import Fraction as Fr
from math import factorial
import flint


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def hook_dim(lam):
    n = sum(lam)
    if n == 0:
        return 1
    conj = [sum(1 for x in lam if x > j) for j in range(lam[0])]
    h = 1
    for i, li in enumerate(lam):
        for j in range(li):
            h *= (li - j - 1) + (conj[j] - i - 1) + 1
    return factorial(n) // h


def perm_sign(images):
    n = len(images)
    seen = [False] * n
    s = 1
    for i in range(n):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = images[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def group_elements(k, groups):
    """all permutations of range(k) (as tuples of images) that permute each list of groups (disjoint) and fix the
    other copies, with their signs"""
    elems = [(tuple(range(k)), 1)]
    for g in groups:
        if len(g) < 2:
            continue
        new = []
        for p in itertools.permutations(g):
            sg = perm_sign(tuple(g.index(x) for x in p))
            for (e, s) in elems:
                e2 = list(e)
                for x, y in zip(g, p):
                    e2[x] = y
                new.append((tuple(e2), s * sg))
        elems = new
    return elems


def act(pi, st):
    a, y = st
    k = len(y)
    a2, y2 = [0] * k, [0] * k
    for t in range(k):
        a2[pi[t]] = a[t]
        y2[pi[t]] = y[t]
    return (tuple(a2), tuple(y2))


def irreps(blocks):
    """all tuples of partitions, one per block"""
    return list(itertools.product(*[list(partitions(len(B))) for B in blocks]))


class Young:
    def __init__(self, k, blocks, lams):
        rows, cols = [], []
        for B, lam in zip(blocks, lams):
            B = sorted(B)
            pos, rws = 0, []
            for li in lam:
                rws.append(B[pos:pos + li])
                pos += li
            rows += rws
            for jj in range(lam[0] if lam else 0):
                col = [rw[jj] for rw in rws if len(rw) > jj]
                if len(col) > 1:
                    cols.append(col)
        self.R = [p for p, s in group_elements(k, rows)]
        self.C = group_elements(k, cols)
        self.dim = 1
        for lam in lams:
            self.dim *= hook_dim(lam)


class Module:
    """the permutation module on a list of states"""

    def __init__(self, states):
        self.states = states
        self.idx = {st: i for i, st in enumerate(states)}
        self._pc = {}

    def img(self, pi, i):
        key = (pi, i)
        v = self._pc.get(key)
        if v is None:
            v = self.idx[act(pi, self.states[i])]
            self._pc[key] = v
        return v

    def apply_e(self, Yg, v):
        """e v for a dict vector v (index -> Fraction)"""
        w = {}
        for (c, sg) in Yg.C:
            for i, x in v.items():
                j = self.img(c, i)
                w[j] = w.get(j, 0) + sg * x
        out = {}
        for r in Yg.R:
            for i, x in w.items():
                if x == 0:
                    continue
                j = self.img(r, i)
                out[j] = out.get(j, 0) + x
        return {i: x for i, x in out.items() if x != 0}

    def check_invariant(self, perms):
        """every permutation maps the state list onto itself"""
        return all(sorted(self.idx[act(p, st)] for st in self.states) == list(range(len(self.states)))
                   for p in perms)


def span_basis(vecs, n):
    """a basis (reduced row echelon rows, as lists of Fractions) of the span of the dict vectors vecs"""
    if not vecs:
        return []
    M = flint.fmpq_mat(len(vecs), n)
    for r, v in enumerate(vecs):
        for i, x in v.items():
            x = Fr(x)
            M[r, i] = flint.fmpq(x.numerator, x.denominator)
    R, rk = M.rref()
    out = []
    for r in range(rk):
        out.append([Fr(int(R[r, c].p), int(R[r, c].q)) for c in range(n)])
    return out
