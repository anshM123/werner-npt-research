"""Test the qubit explanation of the t = -1 local-search zeros at n = 4.
For a zero (u, v): Q_i = support of rho^u_i + rho^v_i (2-dim); write u, v in the product basis of the Q_i; then
  (1) the full-weight Pauli part of P_u - P_v (all 81 sigma_mu, mu in {x,y,z}^4) should vanish;
  (2) test v ~ J^(x)4 conj(u) with J = sigma_y (spin flip), i.e. |<v, Y conj(u)>| = 1.
usage: python tm1_qubit.py file.npy [...]"""
import sys, itertools, numpy as np

SX = np.array([[0, 1], [1, 0]], complex); SY = np.array([[0, -1j], [1j, 0]]); SZ = np.diag([1.0 + 0j, -1.0]); I2 = np.eye(2)
PA = [I2, SX, SY, SZ]

def red_site(z, n, i):
    Z = np.moveaxis(z.reshape((3,) * n), i, 0).reshape(3, -1)
    return Z @ Z.conj().T

if __name__ == "__main__":
    for fn in sys.argv[1:]:
        U = np.load(fn); n = 4; u, v = U[:, 0], U[:, 1]
        Bs = []
        for i in range(n):
            w, V = np.linalg.eigh(red_site(u, n, i) + red_site(v, n, i)); Bs.append(V[:, 1:])      # top-2 eigvecs, 3 x 2
        B = Bs[0]
        for i in range(1, n): B = np.kron(B, Bs[i])                                                    # 81 x 16 isometry
        uq = B.conj().T @ u; vq = B.conj().T @ v
        leak = max(np.linalg.norm(u - B @ uq), np.linalg.norm(v - B @ vq))
        C = np.outer(uq, uq.conj()) - np.outer(vq, vq.conj())
        fw = 0.0; tot = 0.0
        for mu in itertools.product(range(4), repeat=n):
            S = PA[mu[0]]
            for k in mu[1:]: S = np.kron(S, PA[k])
            c = np.trace(S @ C) / 2 ** n
            tot += abs(c) ** 2
            if 0 not in mu: fw += abs(c) ** 2
        Y = SY
        for _ in range(n - 1): Y = np.kron(Y, SY)
        sf = abs(np.vdot(vq, Y @ uq.conj()))
        print(f"{fn}: leak outside Q {leak:.1e}; full-weight Pauli weight of P_u - P_v {fw:.2e} (total {tot:.3f}); "
              f"|<v, Y conj u>| = {sf:.6f}; <u, Y conj u> = {abs(np.vdot(uq, Y @ uq.conj())):.2e}")
