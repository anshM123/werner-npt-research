# Root locations of F_Z(t) = sum_J t^{|J|} ||Tr_J Z||^2 for random rank-2 Z on (C^3)^{(x)n}
import numpy as np, itertools, sys
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 0)
n = int(sys.argv[1]); d = 3
def coeffs_multi(Z):
    T = Z.reshape([d]*n + [d]*n)
    c = {}
    for J in itertools.product([0,1], repeat=n):
        Jl = [i for i in range(n) if J[i]]
        X = T
        # trace over sites in Jl (pair ket axis i with bra axis n+i), do from highest index
        X = Z.reshape([d]*(2*n))
        axes_left = list(range(n))
        cur = X; m = n
        for i in sorted(Jl, reverse=True):
            cur = np.trace(cur, axis1=i, axis2=i+m); m -= 1
        c[J] = np.sum(np.abs(cur)**2)
    return c
nonreal = 0; inside = 0; tot = 0; worst = 1e9
for trial in range(int(sys.argv[3]) if len(sys.argv)>3 else 200):
    kind = trial % 3
    if kind == 0:
        A = rng.normal(size=(d**n,2)) + 1j*rng.normal(size=(d**n,2)); B = rng.normal(size=(d**n,2)) + 1j*rng.normal(size=(d**n,2))
        Z = A @ B.conj().T
    elif kind == 1:
        A = rng.normal(size=(d**n,2)) + 1j*rng.normal(size=(d**n,2)); Q,_ = np.linalg.qr(A); Z = Q @ Q.conj().T
    else:
        A = rng.normal(size=(d**n,2)); Q,_ = np.linalg.qr(A); Z = Q @ Q.T
    c = coeffs_multi(Z)
    m = np.zeros(n+1)
    for J,v in c.items(): m[sum(J)] += v
    r = np.roots(m[::-1])
    tot += 1
    if np.max(np.abs(r.imag)) > 1e-8: nonreal += 1
    # any root with |t| < 1/2 ?
    if np.min(np.abs(r)) < 0.5 - 1e-9: inside += 1
    worst = min(worst, np.min(np.abs(r)))
print(f"n={n}: {tot} trials, nonreal-rooted {nonreal}, root inside |t|<1/2: {inside}, min |root| = {worst:.4f}")
