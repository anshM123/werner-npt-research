"""Standalone exact verifier for the qutrit pair certificate (pure Python + fractions, no floating point).

Claim: at beta = 83/200, the two-site qutrit Werner product rho_beta (x) rho_beta (rho_beta = I - beta F on C^3 (x) C^3;
Alice holds sites a1 a2, Bob holds b1 b2) is a convex combination of the twirled operators
    X_j = (A_j (x) B_j) E_9 (A_j (x) B_j)^T / Tr(...),      E_9 = I - F/2 on C^9 (x) C^9,
with the real rational filters A_j, B_j : C^9 -> C^9 stored in pair_certificate_83_200.json.
Twirled two-site operators are determined by their four sector weights (ss, sa, as, aa), so it suffices to check that
the target sector weights are a convex combination of the generators' sector weights.
Consequence (by the two-copy theorem in dimension 9): rho_{83/200}^{(x)4} is one-copy undistillable, i.e. qutrit Werner
states with alpha <= 83/200 are four-copy undistillable.
usage: verify_pair_certificate.py [pair_certificate_83_200.json]
"""
import itertools
import json
import sys
from fractions import Fraction as Fr

path = sys.argv[1] if len(sys.argv) > 1 else 'pair_certificate_83_200.json'
cert = json.load(open(path))
beta = Fr(*cert['beta'])
den = cert['denominator']
d, N = 3, 9

# basis index helpers: Alice index a = 3 a1 + a2, Bob index b = 3 b1 + b2, joint index (a, b)


def sector_weight(X, sector):
    """Tr[X Pi_sector] for X given as dict ((a, b), (a', b')) -> value; Pi = P^{x1}_{site1} (x) P^{x2}_{site2}.
    P^s = (1 + F)/2, P^a = (1 - F)/2 on (a_i, b_i);  F swaps a_i <-> b_i."""
    s1, s2 = sector
    tot = Fr(0)
    for (row, col), v in X.items():
        (a, b), (ap, bp) = row, col
        a1, a2 = divmod(a, 3); b1, b2 = divmod(b, 3)
        ap1, ap2 = divmod(ap, 3); bp1, bp2 = divmod(bp, 3)
        # <col| Pi |row> component: Pi_{(ap,bp),(a,b)} = prod over sites of P_{(ap_i, bp_i),(a_i, b_i)}
        def P(kind, x, y, xp, yp):        # <xp yp| P |x y>
            ident = Fr(1) if (xp == x and yp == y) else Fr(0)
            swap = Fr(1) if (xp == y and yp == x) else Fr(0)
            return (ident + swap) / 2 if kind == 's' else (ident - swap) / 2
        tot += v * P(s1, a1, b1, ap1, bp1) * P(s2, a2, b2, ap2, bp2)
    return tot


SECTORS = [('s', 's'), ('s', 'a'), ('a', 's'), ('a', 'a')]


def profile(X):
    w = [sector_weight(X, s) for s in SECTORS]
    t = sum(w)
    return [x / t for x in w]


# target: rho_beta (x) rho_beta with rho = I - beta F on each site
T = {}
for a, b, ap, bp in itertools.product(range(9), repeat=4):
    a1, a2 = divmod(a, 3); b1, b2 = divmod(b, 3); ap1, ap2 = divmod(ap, 3); bp1, bp2 = divmod(bp, 3)
    def rho(x, y, xp, yp):
        return (Fr(1) if (x == xp and y == yp) else Fr(0)) - beta * (Fr(1) if (xp == y and yp == x) else Fr(0))
    v = rho(a1, b1, ap1, bp1) * rho(a2, b2, ap2, bp2)
    if v != 0:
        T[((ap, bp), (a, b))] = v
target = profile(T)
print('target profile:', [str(x) for x in target])


def core_image(A, B):
    """(A (x) B) E_9 (A (x) B)^T with E_9 = I - F/2 on C^9 (x) C^9 (as a sparse dict)."""
    A = [[Fr(x, den) for x in row] for row in A]
    B = [[Fr(x, den) for x in row] for row in B]
    # columns of A (x) B indexed by (k, l): vector A[:,k] (x) B[:,l]
    # E_9 = sum_{k,l} |kl><kl| - 1/2 sum_{k,l} |kl><lk|
    X = {}
    cols = {(k, l): [(a, b, A[a][k] * B[b][l]) for a in range(9) for b in range(9) if A[a][k] != 0 and B[b][l] != 0]
            for k in range(9) for l in range(9)}
    for k in range(9):
        for l in range(9):
            for (k2, l2, coef) in ((k, l, Fr(1)), (l, k, Fr(-1, 2))):
                u = cols[(k, l)]
                v = cols[(k2, l2)]
                for (a, b, x) in u:
                    for (ap, bp, y) in v:
                        key = ((a, b), (ap, bp))
                        X[key] = X.get(key, Fr(0)) + coef * x * y
    return X


profs = []
for g in cert['generators']:
    assert g['kind'] == 'core'
    profs.append(profile(core_image(g['A'], g['B'])))
    print('generator profile:', [float(x) for x in profs[-1]])

# solve sum_j lam_j prof_j = target, sum lam_j = 1 exactly (Gaussian elimination over Fractions)
k = len(profs)
rows = [[profs[j][r] for j in range(k)] + [target[r]] for r in range(4)] + [[Fr(1)] * k + [Fr(1)]]
M = [r[:] for r in rows]
piv_cols = []
r0 = 0
for c in range(k):
    p = next((r for r in range(r0, len(M)) if M[r][c] != 0), None)
    if p is None:
        continue
    M[r0], M[p] = M[p], M[r0]
    for r in range(len(M)):
        if r != r0 and M[r][c] != 0:
            f = M[r][c] / M[r0][c]
            M[r] = [x - f * y for x, y in zip(M[r], M[r0])]
    piv_cols.append(c)
    r0 += 1
lam = [Fr(0)] * k
for i, c in enumerate(piv_cols):
    lam[c] = M[i][k] / M[i][c]
residual_ok = all(sum(rows[r][j] * lam[j] for j in range(k)) == rows[r][k] for r in range(len(rows)))
print('weights (float):', [float(x) for x in lam])
print('identity exact:', residual_ok, '| all weights >= 0:', all(x >= 0 for x in lam))
print('CERTIFICATE VALID:', residual_ok and all(x >= 0 for x in lam))
