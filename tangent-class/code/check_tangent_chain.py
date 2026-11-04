import numpy as np
from fock_tangent import *
rng = np.random.default_rng(1)
worst = 0
for n in (2, 3, 4):
    for trial in range(3):
        U, a = rand_tangent(rng, n)
        O = Ops(U, a); D = O.D
        v = rng.normal(size=(2, D)) + 1j * rng.normal(size=(2, D))
        u = tangent_range(U, a, n); X = np.outer(u[0], v[0].conj()) + np.outer(u[1], v[1].conj())
        F = O.F_of_v(v)
        qd = q_direct(X, n); qf = np.vdot(F, O.Q @ F).real
        # gauge check sum Abar + a a^* = I
        g = sum(O.Ab) + np.outer(a, a.conj())
        # chain form: sum_c <F_c|I + Vac - T + 1/2 al^* al|F_c> - sqrt2 Re sum_c <al F_{c+1}, <a|F_c>
        Dc = np.eye(2 * D) + O.Vac - O.T + 0.5 * O.al.conj().T @ O.al
        tot = 0
        Fc = [np.where(np.concatenate([O.exc, O.exc]) == c, F, 0) for c in range(n + 1)]
        for c in range(n + 1): tot += np.vdot(Fc[c], Dc @ Fc[c]).real
        for c in range(n): tot -= np.sqrt(2) * np.real(np.vdot(O.al @ Fc[c + 1], O.abra @ Fc[c]))
        # second form: 2q = <F|4Vac - 2 sum Delta|F> + ||Afr F||^2, Delta_t = beta_t^* beta_t - Abar_t (x) P_{1,t}
        Ex = O.AbU - O.Vac
        val2 = np.vdot(F, (4 * O.Vac - 2 * (O.T - Ex)) @ F).real + np.linalg.norm(O.Afr @ F) ** 2
        e1 = abs(qf - qd) / abs(qd); e2 = abs(tot - qd) / abs(qd); e3 = abs(val2 / 2 - qd) / abs(qd)
        worst = max(worst, e1, e2, e3, np.abs(g - np.eye(2)).max())
        print(n, trial, f"q={qd:.6f} fock={qf:.6f} chain={tot:.6f} delta-form={val2/2:.6f}")
print("worst rel err", worst); assert worst < 1e-10; print("PASS")
