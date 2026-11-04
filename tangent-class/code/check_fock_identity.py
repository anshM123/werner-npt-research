"""Check of the Fock identity for tangent planes (Proposition 1), with arbitrary (unnormalised) range blocks U_s,
vacuum amplitudes a and co-range v, against explicit partial traces.  a = 0 is the single-excitation case.
  q_n(X) = <F| sum_u Abar_u (x) (2|0><0|_u + Pi_{E_u}) |F> - sum_s ||Z_s F||^2 + 1/2 || sum_s Z_s F - sqrt2 (<a| (x) I) F ||^2,
X = sum_r u_r v_r^*, F = sum_r |r> (x) D^{1/2} v_r, D = (x) diag(1/2, 1, 1), Abar_s = (U_s^* U_s)^T,
Z_s F = sum_{r,z} conj(U_s[z,r]) <z|_s f_r with |0> placed at copy s."""
import numpy as np
from common import strings, q_direct, tangent_range, d
rng = np.random.default_rng(2026)


def lower(f, s, n, z):
    T = f.reshape([d] * n); g = np.zeros_like(T)
    src = [slice(None)] * n; src[s] = z; dst = [slice(None)] * n; dst[s] = 0
    g[tuple(dst)] = T[tuple(src)]; return g.reshape(-1)


worst = 0.0
for n in (2, 3, 4):
    st = strings(n); D = d ** n
    Dh = np.array([np.prod([np.sqrt(0.5) if c == 0 else 1.0 for c in s]) for s in st])
    for trial in range(4):
        U = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(n)]
        a = (rng.normal(size=2) + 1j * rng.normal(size=2)) * (0 if trial == 0 else 1)
        u = tangent_range(U, a, n)
        v = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D))
        X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj())
        f = [Dh * v[r] for r in range(2)]
        Ab = [(U[s].conj().T @ U[s]).T for s in range(n)]
        one = 0.0
        for i, s_ in enumerate(st):
            fv = np.array([f[0][i], f[1][i]])
            H = sum(Ab[w] * (2 if s_[w] == 0 else 1) for w in range(n))
            one += (fv.conj() @ H @ fv).real
        Z = [sum(U[s][z - 1, r].conj() * lower(f[r], s, n, z) for r in range(2) for z in (1, 2)) for s in range(n)]
        Fa = a[0].conj() * f[0] + a[1].conj() * f[1]
        val = one - sum(np.sum(abs(z) ** 2) for z in Z) + 0.5 * np.sum(abs(sum(Z) - np.sqrt(2) * Fa) ** 2)
        qd = q_direct(X, n); err = abs(val - qd) / max(1.0, abs(qd)); worst = max(worst, err)
        print(f"n={n} trial {trial} ({'a = 0' if trial == 0 else 'a != 0'}): q_n = {qd:.10f}, Fock form = {val:.10f}, rel. error {err:.1e}")
print(f"RESULT: maximal relative error {worst:.1e}")
assert worst < 1e-10
print("PASS")
