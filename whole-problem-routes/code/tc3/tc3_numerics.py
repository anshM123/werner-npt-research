"""Numerics quoted around Theorem 7.29 (floating point).

1. Two copies: the best constant k in lambda_min(P^{(x)2}(yy^T)) >= (1/16)(1 + k (tau_1 + tau_2)) is 5, attained at the
   maximally entangled vector (grid over the Schmidt simplex with the block formula of Theorem 7.2).
2. Three copies: projected gradient search for the best such constant; the smallest value found is 5 to search
   accuracy, approached at e (x) Omega-hat (which gives exactly 5).
3. Conditioning on one copy loses thickness: y = (1/2)[e1 (e1 e1 + e2 e3) + e2 (e1 e2 + e3 e3)] is thin at copy 1 and
   has tau = 1/4 at copies 2 and 3, but every conditional vector (f^T (x) I)y with f in the plane of copy 1 has Schmidt
   coefficients (1/sqrt2, 1/sqrt2, 0).  Its least output eigenvalue 0.0774 is far above the TC value 0.0352, while the
   entanglement-breaking bound of Proposition 7.8 with the two-copy bound on the conditional vectors only gives 1/32.
Ends with PASS or FAIL."""
import numpy as np
import torch
from scipy.optimize import minimize

torch.set_default_dtype(torch.float64)
a = 0.75
ok = True

# ---------------------------------------------------------------- 1
best = np.inf
for u in np.linspace(0.0, np.pi / 2, 601):
    for v in np.linspace(0.0, np.pi / 2, 601):
        s = np.array([np.sin(u) * np.cos(v), np.sin(u) * np.sin(v), np.cos(u)])
        tau = np.min(s ** 2)
        if tau < 1e-6:
            continue
        MA = np.diag(a ** 2 - (a - 1 / 8) * s ** 2) + np.outer(s, s) / 8
        vals = [np.linalg.eigvalsh(MA)[0]] + [a ** 2 - (a / 2) * (s[i] ** 2 + s[j] ** 2) - s[i] * s[j] / 8
                                              for i in range(3) for j in range(i + 1, 3)]
        best = min(best, (16 * min(vals) - 1) / (2 * tau))
s = np.ones(3) / np.sqrt(3)
MA = np.diag(a ** 2 - (a - 1 / 8) * s ** 2) + np.outer(s, s) / 8
lme = min([np.linalg.eigvalsh(MA)[0]] + [a ** 2 - (a / 2) * (s[i] ** 2 + s[j] ** 2) - s[i] * s[j] / 8
                                         for i in range(3) for j in range(i + 1, 3)])
kme = (16 * lme - 1) / (2 / 3)
print(f'1. two copies: least linear constant on the grid {best:.6f}; at the maximally entangled vector {kme:.12f} '
      f'(lambda_min = {lme:.12f} = 13/48)')
ok &= best >= 5 - 1e-9 and abs(kme - 5) < 1e-12 and abs(lme - 13 / 48) < 1e-14


# ---------------------------------------------------------------- 2
def pmap_t(X, m):
    T = X.reshape([3] * (2 * m))
    for s in range(m):
        T = torch.movedim(T, (s, m + s), (0, 1))
        sh = T.shape
        Mi = T.reshape(3, 3, -1)
        tr = torch.einsum('iir->r', Mi)
        T = torch.movedim((-0.25 * (Mi + Mi.transpose(0, 1)) + a * torch.eye(3).unsqueeze(-1) * tr).reshape(sh), (0, 1), (s, m + s))
    return T.reshape(3 ** m, 3 ** m)


def taus(v):
    return torch.stack([torch.linalg.eigvalsh(torch.movedim(v.reshape(3, 3, 3), s, 0).reshape(3, -1) @
                                              torch.movedim(v.reshape(3, 3, 3), s, 0).reshape(3, -1).T)[0] for s in range(3)])


def obj(z):
    zt = torch.tensor(z, requires_grad=True)
    y = zt / zt.norm()
    W = pmap_t(torch.outer(y, y), 3)
    val = (64 * torch.linalg.eigvalsh((W + W.T) / 2)[0] - 1) / taus(y).clamp(min=1e-9).sum()
    val.backward()
    return val.item(), zt.grad.numpy()


rng = np.random.default_rng(11)
best3 = np.inf
e = np.zeros(27)
e[0] = 1
for st in range(24):
    y0 = [rng.standard_normal(27), e + 0.3 * rng.standard_normal(27),
          np.kron(np.eye(3)[0], np.eye(3).reshape(-1) / np.sqrt(3)) + 0.2 * rng.standard_normal(27)][st % 3]
    best3 = min(best3, minimize(obj, y0, jac=True, method='L-BFGS-B', options={'maxiter': 4000}).fun)
print(f'2. three copies: smallest linear constant found {best3:.5f}')
ok &= 4.99 < best3 < 5.05

# ---------------------------------------------------------------- 3
def pmap(X, m):
    T = X.reshape([3] * (2 * m))
    for s in range(m):
        T = np.moveaxis(T, [s, m + s], [0, 1])
        sh = T.shape
        Mi = T.reshape(3, 3, -1)
        tr = np.einsum('iir->r', Mi)
        T = np.moveaxis((-0.25 * (Mi + Mi.transpose(1, 0, 2)) + a * np.eye(3)[:, :, None] * tr).reshape(sh), [0, 1], [s, m + s])
    return T.reshape(3 ** m, 3 ** m)


E = np.eye(3)
Y1 = np.array([[1, 0, 0], [0, 0, 1], [0, 0, 0.]])
Y2 = np.array([[0, 1, 0], [0, 0, 0], [0, 0, 1.]])
y = (np.kron(E[0], Y1.reshape(-1)) + np.kron(E[1], Y2.reshape(-1))) / 2
tv = [np.linalg.eigvalsh(np.moveaxis(y.reshape(3, 3, 3), s, 0).reshape(3, -1) @ np.moveaxis(y.reshape(3, 3, 3), s, 0).reshape(3, -1).T)[0] for s in range(3)]
lam = np.linalg.eigvalsh(pmap(np.outer(y, y), 3))[0]
tc = np.prod([0.25 + x / 2 for x in tv])
sch = []
for th in np.linspace(0, np.pi, 13):
    f = np.cos(th) * E[0] + np.sin(th) * E[1]
    cond = np.tensordot(f, y.reshape(3, 3, 3), axes=([0], [0]))
    sv = np.linalg.svd(cond / np.linalg.norm(cond), compute_uv=False)
    sch.append(np.abs(sv - np.array([1, 1, 0]) / np.sqrt(2)).max())
print(f'3. pencil input: tau = {np.round(tv, 6)}, lambda_min = {lam:.6f}, TC value {tc:.6f}, '
      f'conditional Schmidt coefficients (1/sqrt2, 1/sqrt2, 0) to {max(sch):.1e}')
ok &= abs(tv[0]) < 1e-12 and abs(tv[1] - 0.25) < 1e-12 and abs(tv[2] - 0.25) < 1e-12 and lam > tc + 0.04 and max(sch) < 1e-12
print('PASS' if ok else 'FAIL')
