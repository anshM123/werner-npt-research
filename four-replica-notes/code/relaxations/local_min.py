"""Local search for min over rank-r projectors P on (C^3)^(x)n of q_n^{(w)}(P) = sum_J (-w)^{|J|} ||Tr_J P||_F^2
(upper bounds for the true minimum, to compare with the relaxation).  usage: python local_min.py n w [rank] [restarts]"""
import sys, itertools, numpy as np
from scipy.optimize import minimize

def make_q(n, w):
    subsets = [J for k in range(n + 1) for J in itertools.combinations(range(n), k)]
    def q(P):
        Pt = P.reshape((3,) * (2 * n)); tot = 0.0
        for J in subsets:
            rows = [chr(97 + s) for s in range(n)]; cols = [chr(97 + n + s) for s in range(n)]
            for s in J: cols[s] = rows[s]
            out = ''.join(rows[s] for s in range(n) if s not in J) + ''.join(cols[s] for s in range(n) if s not in J)
            R = np.einsum(''.join(rows) + ''.join(cols) + '->' + out, Pt)
            tot += (-w) ** len(J) * np.sum(np.abs(R) ** 2)
        return float(tot)
    return q

def run(n, w, rank=2, restarts=20, seed=0, spec=None):
    q = make_q(n, w); D = 3 ** n
    S = np.diag(spec) if spec is not None else np.eye(rank)
    rng = np.random.default_rng(seed)
    def f(x):
        U = (x[:D * rank] + 1j * x[D * rank:]).reshape(D, rank)
        Q, _ = np.linalg.qr(U)
        return q(Q @ S @ Q.conj().T)
    best = np.inf
    for r in range(restarts):
        x0 = rng.standard_normal(2 * D * rank)
        res = minimize(f, x0, method='L-BFGS-B', options=dict(maxiter=3000))
        best = min(best, res.fun)
    return best

if __name__ == "__main__":
    n = int(sys.argv[1]); w = float(sys.argv[2]); rank = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    restarts = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    spec = [float(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else None
    b = run(n, w, rank, restarts, spec=spec)
    if spec is not None:
        print(f"n={n} w={w} spec={spec}: local-search min {b:.10f};  product value {(sum(x*x for x in spec) - w * sum(spec) ** 2) * (1 - w) ** (n - 1):.10f}")
        sys.exit(0)
    prod = (rank - w * rank ** 2) * (1 - w) ** (n - 1)
    print(f"n={n} w={w} rank={rank}: local-search min {b:.10f};  product P_1 (x) |0><0|^(n-1): {prod:.10f}")
