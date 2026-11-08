"""LP relaxation of the sharp rank-two bound via purity profiles. C >= 0 of rank <= 2 on n qudits (d = 3), Tr C = 1,
f(J) = Tr (C_J)^2 (marginal purities), purified by a qubit R. Constraints:
  f(empty) = 1, f([n]) = P2 in [1/2, 1];
  Rains shadows of C:            Sh(S)   = sum_{J subset S} (-1)^|J| f(J) >= 0;
  shadows including the qubit:   ShR(S)  = sum_{J subset S} (-1)^|J| [f(J) - f(n minus J)] >= 0;
  qubit purification:             f(J) >= f(n minus J)/2 and f(n minus J) >= f(J)/2;   d^-|J| <= f(J) <= 1.
Minimise 2^n q_n - (2 P2 - 1) = sum_J (-1)^|J| 2^(n-|J|) f(n minus J) - (2 P2 - 1).  >= 0 would prove the sharp bound
from these constraints."""
import itertools, sys
import numpy as np
from scipy.optimize import linprog
d = 3
for n in (3, 4, 5, 6, 7):
    subs = [frozenset(S) for k in range(n + 1) for S in itertools.combinations(range(n), k)]
    idx = {S: i for i, S in enumerate(subs)}; full = frozenset(range(n)); m = len(subs)
    worst = np.inf
    for P2 in np.linspace(0.5, 1.0, 11):
        A_ub, b_ub = [], []
        for S in subs:
            row = np.zeros(m)
            for J in subs:
                if J <= S:
                    row[idx[J]] += (-1) ** len(J)
            if S:
                A_ub.append(-row); b_ub.append(0.0)               # Sh(S) >= 0
            rowR = np.zeros(m)
            for J in subs:
                if J <= S:
                    rowR[idx[J]] += (-1) ** len(J); rowR[idx[full - J]] -= (-1) ** len(J)
            A_ub.append(-rowR); b_ub.append(0.0)                  # ShR(S) >= 0
            r1 = np.zeros(m); r1[idx[full - S]] -= 1; r1[idx[S]] += 0.5
            A_ub.append(r1); b_ub.append(0.0)                     # f(S^c) >= f(S)/2
        c = np.zeros(m)
        for J in subs:
            c[idx[full - J]] += (-1) ** len(J) * 2 ** (n - len(J))
        bounds = [(d ** (-len(S)), 1.0) for S in subs]
        A_eq = np.zeros((2, m)); A_eq[0, idx[frozenset()]] = 1; A_eq[1, idx[full]] = 1
        res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=A_eq, b_eq=[1.0, P2], bounds=bounds, method='highs')
        if res.status == 0:
            worst = min(worst, res.fun - (2 * P2 - 1))
    print(f'n = {n}: LP minimum of 2^n q_n - (2 Tr C^2 - (Tr C)^2) over the relaxation: {worst:+.6f}', flush=True)
