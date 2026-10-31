"""Section 1, remark "No basis average works at zeros".

For the zero C = |0^n><0^n| - (-1)^n |1^n><1^n| of q_n(1/2, .), the trace identity
    sum_k Delta_s(e_k) = 2 q_n(C) - (1/2) sum_k q_{n-1}(<e_k|C|e_k>)
is at most -(1/3) 2^{-(n-2)} for every copy s and every orthonormal basis (e_k) of that copy, because
q_{n-1}(<e_k|C|e_k>) = 2^{-(n-1)} (|<e_k|0>|^2 + |<e_k|1>|^2)^2.  The script checks the block formula and the bound on
random bases for n = 2..5, and that the compression w = |2> gives Delta = 0 (so QC itself holds at C).
"""
import numpy as np
from werner import q, contract

d = 3
rng = np.random.default_rng(0)
for n in range(2, 6):
    e0 = np.zeros(d ** n)
    e0[0] = 1
    e1 = np.zeros(d ** n)
    e1[sum(d ** k for k in range(n))] = 1
    C = np.outer(e0, e0) - (-1) ** n * np.outer(e1, e1)
    qn = q(C, n, 0.5)
    worst, formula_err = -np.inf, 0.0
    for trial in range(200):
        U, _ = np.linalg.qr(rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d)))
        for s in range(n):
            blocks = [q(contract(C, n, s, U[:, k], U[:, k]), n - 1, 0.5) for k in range(d)]
            pred = [2.0 ** (-(n - 1)) * (abs(U[0, k]) ** 2 + abs(U[1, k]) ** 2) ** 2 for k in range(d)]
            formula_err = max(formula_err, max(abs(x - y) for x, y in zip(blocks, pred)))
            worst = max(worst, 2 * qn - 0.5 * sum(blocks))
    # QC at C: compress the unused direction |2> at copy 0
    w = np.array([0, 0, 1.0])
    P = np.eye(d) - np.outer(w, w)
    Pn = np.kron(P, np.eye(d ** (n - 1)))
    delta_w = qn - q(Pn @ C @ Pn, n, 0.5)
    print(f"n = {n}: q_n(C) = {qn:+.1e}; block formula error {formula_err:.1e}; "
          f"max over bases and copies of sum_k Delta_s(e_k) = {worst:+.5f} <= -(1/3)2^-(n-2) = {-(2.0 ** -(n - 2)) / 3:+.5f}; "
          f"Delta for w = |2>: {delta_w:+.1e}")
