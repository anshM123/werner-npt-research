"""Exact constructive check of the heterogeneous cone lemma (pipage-rounding proof).

Claim P(n, m): for s in [0, inf)^n with sum_i (1 - s_i)_+ <= m,
    Phi_s(z) = prod_i (s_i + z_i)  is a nonnegative combination of generators  z^E (1+z)^C  (E, C disjoint, |E| <= m).
The proof is constructive; this script runs the construction in exact rational arithmetic, returns the weights, and
verifies (i) all weights >= 0, (ii) every generator has |E| <= m, (iii) the polynomial identity holds exactly.
Recursion (sites outside the current list get factor 1):
  Step 0: s_i > 1           ->  Phi = (s_i - 1) Phi[without i] + Phi[s_i := 1]
  base  : len <= m          ->  expand s_i + z_i = s_i (1 + z_i) + (1 - s_i) z_i
  Step 1: budget < m        ->  lower one s_i towards 0: Phi = (s_i - s_i') Phi[without i] + Phi[s_i := s_i']
  Step 2: budget == m, two fractional i, j -> Phi = lam Phi(s(tau-)) + (1-lam) Phi(s(tau+)) + c Phi[without i, j]
"""
import itertools
import random
import sys
from fractions import Fraction as Fr


def construct(sites, s, m, memo):
    """sites: tuple of site labels; s: dict label -> Fraction.  Returns dict (E, C) -> weight (frozensets)."""
    key = (sites, tuple(s[i] for i in sites), m)
    if key in memo:
        return memo[key]
    out = {}

    def acc(dct, c):
        for g, w in dct.items():
            out[g] = out.get(g, 0) + c * w
    sites = tuple(sites)
    # Step 0
    big = [i for i in sites if s[i] > 1]
    if big:
        i = big[0]
        rest = tuple(x for x in sites if x != i)
        acc(construct(rest, s, m, memo), s[i] - 1)
        s2 = dict(s); s2[i] = Fr(1)
        acc(construct(sites, s2, m, memo), Fr(1))
        memo[key] = out
        return out
    if len(sites) <= m:
        # expand directly
        for choice in itertools.product((0, 1), repeat=len(sites)):
            w = Fr(1); E = []; C = []
            for i, ch in zip(sites, choice):
                if ch == 0:
                    w *= s[i]; C.append(i)
                else:
                    w *= (1 - s[i]); E.append(i)
            if w != 0:
                g = (frozenset(E), frozenset(C))
                out[g] = out.get(g, 0) + w
        memo[key] = out
        return out
    budget = sum(1 - s[i] for i in sites)
    if budget < m:
        # lower the first site with s_i > 0 as far as needed (or to 0)
        i = next(x for x in sites if s[x] > 0)
        new_si = max(Fr(0), s[i] - (m - budget))
        rest = tuple(x for x in sites if x != i)
        acc(construct(rest, s, m, memo), s[i] - new_si)
        s2 = dict(s); s2[i] = new_si
        acc(construct(sites, s2, m, memo), Fr(1))
        memo[key] = out
        return out
    assert budget == m, (budget, m)
    frac = [i for i in sites if 0 < s[i] < 1]
    if not frac:
        E = frozenset(i for i in sites if s[i] == 0)
        C = frozenset(i for i in sites if s[i] == 1)
        assert len(E) == m
        out[(E, C)] = Fr(1)
        memo[key] = out
        return out
    assert len(frac) >= 2
    i, j = frac[0], frac[1]
    tm = -min(s[i], 1 - s[j])
    tp = min(1 - s[i], s[j])
    lam = tp / (tp - tm)
    sm = dict(s); sm[i] = s[i] + tm; sm[j] = s[j] - tm
    sp = dict(s); sp[i] = s[i] + tp; sp[j] = s[j] - tp
    rest = tuple(x for x in sites if x not in (i, j))
    acc(construct(sites, sm, m, memo), lam)
    acc(construct(sites, sp, m, memo), 1 - lam)
    acc(construct(rest, s, m, memo), lam * tm * tm + (1 - lam) * tp * tp)
    memo[key] = out
    return out


def poly_of_generators(weights, n):
    P = {}
    for (E, C), w in weights.items():
        Cl = sorted(C)
        for k in range(len(Cl) + 1):
            for sub in itertools.combinations(Cl, k):
                S = frozenset(E) | frozenset(sub)
                P[S] = P.get(S, 0) + w
    return P


def target_poly(s, n):
    P = {}
    for k in range(n + 1):
        for S in itertools.combinations(range(n), k):
            w = Fr(1)
            for i in range(n):
                if i not in S:
                    w *= s[i]
            P[frozenset(S)] = w
    return P


random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
checked = 0
for trial in range(400):
    n = random.randint(1, 7)
    m = random.randint(0, n)
    # random s with budget exactly m or below, entries rational; include some s_i > 1 and some 0/1
    s = {}
    for i in range(n):
        r = random.random()
        if r < 0.1:
            s[i] = Fr(random.randint(11, 30), 10)
        elif r < 0.2:
            s[i] = Fr(random.choice([0, 1]))
        else:
            s[i] = Fr(random.randint(1, 99), 100)
    budget = sum(max(Fr(0), 1 - s[i]) for i in range(n))
    if budget > m:
        # scale the deficits down to fit the budget exactly
        scale = Fr(m) / budget if budget > 0 else Fr(0)
        for i in range(n):
            if s[i] < 1:
                s[i] = 1 - (1 - s[i]) * scale
    W = construct(tuple(range(n)), s, m, {})
    assert all(w >= 0 for w in W.values()), 'negative weight'
    assert all(len(E) <= m and not (E & C) for (E, C) in W), 'bad generator'
    P1 = poly_of_generators({g: w for g, w in W.items() if w != 0}, n)
    P2 = target_poly(s, n)
    keys = set(P1) | set(P2)
    assert all(P1.get(k, 0) == P2.get(k, 0) for k in keys), 'identity fails'
    checked += 1
print(f'exact constructive check passed on {checked} random instances (n <= 7, all m, budgets <= m, '
      f'including s_i > 1, s_i in {{0,1}})')
