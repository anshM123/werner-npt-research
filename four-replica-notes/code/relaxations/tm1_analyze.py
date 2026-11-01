"""Structure of local-search zeros (u, v) of q_n(P_u + t P_v).
For each zero: q, the kernel values g(u,u), g(v,v), g(u,v) (g(x,y) = <x(x)y| (x)_i (2F_i - I) |x(x)y>), per-site reduced
spectra of u and v, per-site overlaps Tr(rho^u_i rho^v_i), local support of span(u, v) at each site, Schmidt spectra across
every cut {i} | rest and {i,j} | rest, and the best rank-one fit of u and v at each site.
usage: python tm1_analyze.py file.npy [...]"""
import sys, itertools, numpy as np

def qform(C, n):
    T = C.reshape((3,) * (2 * n)); tot = 0.0
    for k in range(n + 1):
        for J in itertools.combinations(range(n), k):
            rows = [chr(97 + s) for s in range(n)]; cols = [chr(97 + n + s) for s in range(n)]
            for s in J: cols[s] = rows[s]
            keep = [s for s in range(n) if s not in J]
            sub = ''.join(rows[s] for s in keep) + ''.join(cols[s] for s in keep)
            R = np.einsum(''.join(rows) + ''.join(cols) + '->' + sub, T)
            tot += (-0.5) ** k * float(np.real(np.vdot(R, R)))
    return tot

def gker(x, y, n):
    """<x(x)y| (x)_i (2F_i - I) |x(x)y> via sum over subsets: 2F - I = sum; use purity form"""
    # g(x, y) = sum_K (-1)^{n-|K|} 2^{|K|} <x(x)y| F_K |x(x)y> = sum_K (-1)^{n-|K|} 2^{|K|} Tr(rho^x_K rho^y_K)
    X = x.reshape((3,) * n); Y = y.reshape((3,) * n); tot = 0.0
    for k in range(n + 1):
        for K in itertools.combinations(range(n), k):
            rest = [s for s in range(n) if s not in K]
            def red(Z):
                M = np.moveaxis(Z, list(K) + rest, list(range(n))).reshape(3 ** len(K), -1)
                return M @ M.conj().T
            val = np.real(np.trace(red(X) @ red(Y))) if k > 0 else 1.0
            tot += (-1) ** (n - k) * 2 ** k * val
    return tot

def red_site(z, n, sites):
    Z = z.reshape((3,) * n); rest = [s for s in range(n) if s not in sites]
    M = np.moveaxis(Z, list(sites) + rest, list(range(n))).reshape(3 ** len(sites), -1)
    return M @ M.conj().T

def schmidt(z, n, sites):
    Z = z.reshape((3,) * n); rest = [s for s in range(n) if s not in sites]
    M = np.moveaxis(Z, list(sites) + rest, list(range(n))).reshape(3 ** len(sites), -1)
    return np.linalg.svd(M, compute_uv=False)

def fmt(a, k=4):
    return '[' + ', '.join(f'{v:.{k}f}' for v in a) + ']'

if __name__ == "__main__":
    for fn in sys.argv[1:]:
        U = np.load(fn); D = U.shape[0]; n = int(round(np.log(D) / np.log(3)))
        u, v = U[:, 0], U[:, 1]
        print(f"== {fn}: n = {n}, <u,v> = {abs(np.vdot(u, v)):.1e}")
        for t in (-1.0, 1.0):
            C = np.outer(u, u.conj()) + t * np.outer(v, v.conj())
            print(f"   q_n(P_u {t:+.0f} P_v) = {qform(C, n):+.3e}")
        guu, gvv, guv = gker(u, u, n), gker(v, v, n), gker(u, v, n)
        print(f"   g(u,u) = {guu:.6f}, g(v,v) = {gvv:.6f}, g(u,v) = {guv:+.6f}")
        for i in range(n):
            ru = red_site(u, n, [i]); rv = red_site(v, n, [i])
            eu = np.linalg.eigvalsh(ru)[::-1]; ev = np.linalg.eigvalsh(rv)[::-1]
            es = np.linalg.eigvalsh(ru + rv)[::-1]
            print(f"   site {i}: spec u {fmt(eu)}, spec v {fmt(ev)}, Tr(rho_u rho_v) {np.real(np.trace(ru @ rv)):.2e}, "
                  f"spec(rho_u + rho_v) {fmt(es)}")
        for k in (1, 2):
            for K in itertools.combinations(range(n), k):
                if k == 2 and n < 4: continue
                print(f"   cut {K}: Schmidt u {fmt(schmidt(u, n, list(K))[:4], 3)}, v {fmt(schmidt(v, n, list(K))[:4], 3)}")
