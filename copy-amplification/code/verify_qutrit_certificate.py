"""Exact check of the qutrit four-copy certificate

    X(8/3)^{(x)4} = Avg[ (1/6561) O^4 + (32/6561) O^3 E + (384/6561) O^2 E^2 + (6144/6561) S E^3 ],

with X(r) = P^s + r P^a on C^3 (x) C^3, O = X(0), E = X(3), S = X(2) (the separable boundary state at d = 3), and Avg
the average over placements of the factors on the four sites.  All operators involved are diagonal in the 16
projectors Pi_A = (x)_{i in A} P^a_i (x)_{i not in A} P^s_i, so the identity holds iff the eigenvalues agree on every
sector A.  Each factor X(r) contributes r on an antisymmetric site and 1 on a symmetric site.
"""
import itertools
from fractions import Fraction as Fr

n = 4
O, E, S = Fr(0), Fr(3), Fr(2)
terms = [(Fr(1, 6561), [O, O, O, O]), (Fr(32, 6561), [O, O, O, E]), (Fr(384, 6561), [O, O, E, E]),
         (Fr(6144, 6561), [S, E, E, E])]


def eig(factors, A):
    v = Fr(1)
    for i, r in enumerate(factors):
        v *= r if i in A else 1
    return v


ok = True
for k in range(n + 1):
    for A in itertools.combinations(range(n), k):
        lhs = Fr(8, 3) ** len(A)
        rhs = Fr(0)
        for w, fac in terms:
            placements = set(itertools.permutations(fac))
            rhs += w * sum(eig(p, A) for p in placements) / len(placements)
        ok &= (lhs == rhs)
        print(f'sector A = {A}: lhs = {lhs}, rhs = {rhs}')
print('identity holds on all 16 sectors:', ok)
print('alpha = (rho - 1)/(rho + 1) at rho = 8/3:', (Fr(8, 3) - 1) / (Fr(8, 3) + 1))
