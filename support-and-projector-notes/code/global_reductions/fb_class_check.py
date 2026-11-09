"""Dense check of the merged sums used in fb_class.py at N = 3, 4 (sites 0, 1 merged)."""
import itertools
import numpy as np
rng = np.random.default_rng(5); om = np.exp(2j * np.pi / 3)
for N in (3, 4):
    for t in (0.4, 0.45):
        lam = 1 - 3 * t; D = 3 ** N
        X = np.array(list(itertools.product(range(3), repeat=N)), dtype=np.int64)
        M = np.triu(rng.integers(0, 3, size=(N, N)), 1); G = (M + M.T) % 3
        ph = sum(G[i, j] * X[:, i] * X[:, j] for i in range(N) for j in range(i + 1, N))
        g = om ** (ph % 3) / np.sqrt(D)
        d = rng.integers(0, 3, N); v = rng.integers(0, 3, N); v[0] = 1
        Zv = lambda u: om ** ((X @ u) % 3)
        C1 = np.outer(Zv(d) * g, g.conj()); C2 = np.outer(Zv((d - v) % 3) * g, (Zv(-v % 3) * g).conj())
        def form(A, B, groups):
            # groups: list of lists of copies forming the sites; weight -t per traced site
            TA = A.reshape((3,) * (2 * N)); TB = B.reshape((3,) * (2 * N)); tot = 0
            for K in itertools.product([0, 1], repeat=len(groups)):
                tr = sorted([c for gi, kk in zip(groups, K) if kk for c in gi], reverse=True)
                a, b, cur = TA, TB, N
                for c in tr:
                    a = np.trace(a, axis1=c, axis2=c + cur); b = np.trace(b, axis1=c, axis2=c + cur); cur -= 1
                tot += (-t) ** sum(K) * np.vdot(a, b)
            return tot
        groups_m = [[0, 1]] + [[s] for s in range(2, N)]
        Qm, Xm = form(C1, C1, groups_m).real, form(C1, C2, groups_m)
        # sums as in fb_class.py
        GA = (X @ G.T) % 3; idn = ~((X != 0) | ((GA + d) % 3 != 0))
        fm = np.where(idn[:, 0] & idn[:, 1], 1 - 9 * t, 1.0) * lam ** idn[:, 2:].sum(1)
        Sm = np.fft.fftn(fm.reshape((3,) * N)).ravel() / 3 ** N
        vi = int(np.ravel_multi_index(tuple(v), (3,) * N)); vneg = int(np.ravel_multi_index(tuple((-v) % 3), (3,) * N))
        print(f"N={N} t={t}: Q' dense {Qm:.6f} vs sum {Sm[0].real:.6f};  |X'| dense {abs(Xm):.6f} vs sums {abs(Sm[vi]):.6f}, {abs(Sm[vneg]):.6f}")
