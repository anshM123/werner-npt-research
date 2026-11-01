"""Local search for min q_n(1/2, C) over Hermitian C = U diag(1, t) U^dagger, U an isometry C^2 -> (C^3)^(x)n.
q(C) = Tr[C G(C)] with G(C) = sum_J (-1/2)^{|J|} (Tr_J C) (x) I_J; Riemannian gradient descent on the Stiefel manifold
(Euclidean gradient 4 G(C) U diag(1,t), QR retraction, Barzilai-Borwein steps with Armijo backtracking), random starts.
usage: python local_herm.py n t [restarts] [iters] [seed]"""
import sys, itertools, numpy as np

def make_G(n):
    subsets = [J for k in range(n + 1) for J in itertools.combinations(range(n), k)]
    D = 3 ** n
    def G(C):
        T = C.reshape((3,) * (2 * n)); out = np.zeros((3,) * (2 * n), complex)
        for J in subsets:
            rows = [chr(97 + s) for s in range(n)]; cols = [chr(97 + n + s) for s in range(n)]
            for s in J: cols[s] = rows[s]
            keep = [s for s in range(n) if s not in J]
            sub_out = ''.join(rows[s] for s in keep) + ''.join(cols[s] for s in keep)
            R = np.einsum(''.join(rows) + ''.join(cols) + '->' + sub_out, T)
            if J:
                # (Tr_J C) (x) I_J back in site order
                ops = [R] + [np.eye(3)] * len(J)
                subs = [sub_out] + [chr(97 + s) + chr(97 + n + s) for s in J]
                full = ''.join(chr(97 + s) for s in range(n)) + ''.join(chr(97 + n + s) for s in range(n))
                R = np.einsum(','.join(subs) + '->' + full, *ops)
            out = out + (-0.5) ** len(J) * R
        return out.reshape(D, D)
    return G

def run(n, t, restarts=10, iters=3000, seed=0):
    G = make_G(n); D = 3 ** n; lam = np.diag([1.0, t])
    rng = np.random.default_rng(seed)
    def f_and_grad(U):
        C = U @ lam @ U.conj().T
        GC = G(C)
        q = float(np.real(np.trace(C @ GC)))
        Z = 4 * GC @ U @ lam
        S = U.conj().T @ Z; grad = Z - U @ ((S + S.conj().T) / 2)
        return q, grad
    def retract(U):
        Q, R = np.linalg.qr(U); return Q * np.sign(np.real(np.diag(R)))
    best = np.inf; bestU = None
    for r in range(restarts):
        U = retract(rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))
        q, g = f_and_grad(U); step = 0.05; prevU = None; prevg = None
        for it in range(iters):
            gn = np.real(np.vdot(g, g))
            if gn < 1e-24: break
            if prevU is not None:
                s = (U - prevU).reshape(-1); y = (g - prevg).reshape(-1)
                sy = np.real(np.vdot(s, y))
                if sy > 1e-30: step = min(max(np.real(np.vdot(s, s)) / sy, 1e-6), 10.0)
            while True:
                Un = retract(U - step * g); qn, gnew = f_and_grad(Un)
                if qn <= q - 1e-4 * step * gn or step < 1e-12: break
                step *= 0.5
            prevU, prevg = U, g
            U, q, g = Un, qn, gnew
        if q < best: best, bestU = q, U
    return best, bestU

if __name__ == "__main__":
    n = int(sys.argv[1]); t = float(sys.argv[2])
    restarts = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    iters = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    best, U = run(n, t, restarts, iters, seed)
    prod = (1 + t * t - 0.5 * (1 + t) ** 2) * 0.5 ** (n - 1)
    print(f"n={n} t={t:+.2f}: local-search min {best:.12f}  (product C_1 (x) |0><0|^(n-1): {prod:.12f}; restarts {restarts})", flush=True)
    np.save(f"local_herm_n{n}_t{t:+.2f}.npy", U)
