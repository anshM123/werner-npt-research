import numpy as np
from ex3 import schur_soft
from ex1 import reduced, kernel_basis
np.set_printoptions(linewidth=220, precision=4, suppress=True)
for k in (4,5):
    for big in (3, 10, 30, 100):
        w = np.ones(k); w[0] = big
        m, x, Sh, lamR = schur_soft(w)
        b = k-2
        blocks = [np.linalg.eigvalsh(Sh[h*b:(h+1)*b, h*b:(h+1)*b]) for h in range(k)]
        m2, Y, F = reduced(list(w))
        K = kernel_basis(m2); Qk,_ = np.linalg.qr(K); P = np.eye(len(F)) - Qk@Qk.T
        ev = np.linalg.eigvalsh(P@F@P + 100*Qk@Qk.T)
        print(f"k={k} w0={big}: x={np.round(x,5)} min eig F {ev[0]:.3e} (/x_min {ev[0]/x.min():.3f}, /(1-x0) {ev[0]/(1-x[0]):.3f}); lam(Sig_hat) {np.linalg.eigvalsh(Sh)[:3]}; blocks h=0: {blocks[0]}, h=1: {blocks[1]}; lamR {lamR:.3f}")
