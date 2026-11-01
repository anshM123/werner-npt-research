"""Plane-class corollary of Theorem 3 (FOUR_COPY_THEOREMS.md, Section 2).

For C = Y M Y* with Y an isometry from C^2 onto a plane V of (C^3)^(x)4 and M any 2 x 2 matrix:
  (a) q_4(C) = q_4(Re C) + q_4(Im C), with Re C = (C + C*)/2 and Im C = (C - C*)/(2i);
  (b) Re C and Im C are Hermitian of rank at most two.
So Theorem 3 gives q_4(C) >= 0.  This script checks (a) and (b) on random non-normal M.
usage: python check_plane_class.py [trials] [seed]"""
import sys
import numpy as np

def q_form(C, n, t=0.5):
    """q_n(C) = sum_J (-t)^|J| ||Tr_J C||_2^2, computed as <C, R(C)> with R = prod_s (id - t Tr_s(.) (x) I_s)"""
    T = C.reshape((3,) * (2 * n)).copy()
    for s in range(n):
        tr = np.trace(T, axis1=s, axis2=n + s)
        shape = [3 if k in (s, n + s) else 1 for k in range(2 * n)]
        T = T - t * np.expand_dims(np.expand_dims(tr, s), n + s) * np.eye(3).reshape(shape)
    return float(np.real(np.vdot(C, T.reshape(C.shape))))

def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    n, D = 4, 81
    rng = np.random.default_rng(seed)
    dev = 0.0; rank = 0; qmin = np.inf
    for _ in range(trials):
        Y = np.linalg.qr(rng.standard_normal((D, 2)) + 1j * rng.standard_normal((D, 2)))[0]
        M = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
        C = Y @ M @ Y.conj().T
        C /= np.linalg.norm(C)
        H = (C + C.conj().T) / 2
        K = (C - C.conj().T) / 2j
        dev = max(dev, abs(q_form(C, n) - q_form(H, n) - q_form(K, n)))
        rank = max(rank, np.linalg.matrix_rank(H, 1e-10), np.linalg.matrix_rank(K, 1e-10))
        qmin = min(qmin, q_form(C, n))
    print(f"{trials} random non-normal M at n = 4: max |q(C) - q(Re C) - q(Im C)| = {dev:.1e}, "
          f"max rank of Re C, Im C = {rank}, min q_4(C)/||C||^2 = {qmin:.4f}")
    ok = dev < 1e-12 and rank <= 2
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
