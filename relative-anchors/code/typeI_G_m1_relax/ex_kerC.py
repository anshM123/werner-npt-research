"""Exploration (exact): the kernel of C (31 x 31) at rational weights.  Is it the same for two weight vectors with the
same (sigma, S)?  If so, the 31 functionals can be reduced to rank(C) combinations by a change of basis that depends
only on (sigma, S)."""
import flint
from reduced_c import ReducedC

Q = flint.fmpq
R = ReducedC()


def kernel(wt, wU):
    num = R.numeric(wt, wU)
    J = len(num['C'])
    Cm = flint.fmpq_mat(J, J, [num['C'][i][j] for i in range(J) for j in range(J)])
    rr, rank = Cm.rref()
    # nullspace from rref
    piv = []
    r = 0
    for col in range(J):
        if r < rank and rr[r, col] != 0:
            piv.append(col)
            r += 1
    free = [cidx for cidx in range(J) if cidx not in piv]
    basis = []
    for f in free:
        v = [Q(0)] * J
        v[f] = Q(1)
        for i, pcol in enumerate(piv):
            v[pcol] = -rr[i, f]
        basis.append(v)
    return rank, piv, free, basis


# same sigma and S, different distributions (u = 4): (3/5)^2 + (4/5)^2 = 1 = (5/13)^2 + (12/13)^2
wt = Q(1)
A = [Q(3, 5), Q(4, 5), Q(1), Q(1)]
B = [Q(1), Q(1), Q(5, 13), Q(12, 13)]
C2 = [Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2), Q(1, 2)]
print('S_A =', sum((x * x for x in A), Q(0)), ' S_B =', sum((x * x for x in B), Q(0)), ' S_C =', sum((x * x for x in C2), Q(0)))
rA, pA, fA, kA = kernel(wt, A)
numB = R.numeric(wt, B)
J = len(numB['C'])
CB = flint.fmpq_mat(J, J, [numB['C'][i][j] for i in range(J) for j in range(J)])
numC = R.numeric(wt, C2)
CC = flint.fmpq_mat(J, J, [numC['C'][i][j] for i in range(J) for j in range(J)])
print('rank C at A:', rA, '; kernel dimension', len(kA))
inB = all(all((CB * flint.fmpq_mat(J, 1, v))[i, 0] == 0 for i in range(J)) for v in kA)
inC = all(all((CC * flint.fmpq_mat(J, 1, v))[i, 0] == 0 for i in range(J)) for v in kA)
print('kernel of C at A is contained in the kernel at B (same sigma, S, u):', inB)
print('... and at C2 (same sigma, S, u = 12, equal weights):', inC)
names = R.gens
print('dependent functionals (free columns) at A:', [names[i] for i in fA])
