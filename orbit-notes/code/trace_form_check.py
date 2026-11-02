"""Control: the orbit property is NOT a property of every unitarily invariant form.
For q(X) = ||X||^2 - c |Tr X|^2 on (C^3)^{(x)2} and a positive semidefinite rank-two C that is thick at both copies,
R = 1 - c |Tr X|^2 / ||X||^2 is minimised on the squeeze orbit where G C G is proportional to a rank-two projector
(|Tr X|^2 = 2 ||X||^2).  Those points are interior points of the orbit.
usage: python trace_form_check.py starts seed c"""
import sys, numpy as np, torch
from scipy.optimize import minimize
torch.set_default_dtype(torch.float64)
STARTS, SEED, c = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]); n = 2; D = 9
torch.manual_seed(SEED)
def R(X): return 1 - c * (abs(torch.trace(X)) ** 2 / torch.vdot(X.reshape(-1), X.reshape(-1)).real)
def herm(x):
    M = torch.zeros(3, 3, dtype=torch.complex128) + torch.diag(x[:3].to(torch.complex128)); iu = torch.triu_indices(3, 3, 1)
    M[iu[0], iu[1]] = torch.complex(x[3:6], x[6:9]); M[iu[1], iu[0]] = torch.complex(x[3:6], -x[6:9]); return M
def squeeze(C, x):
    G = torch.kron(torch.linalg.matrix_exp(herm(x[:9])), torch.linalg.matrix_exp(herm(x[9:])))
    return G @ C @ G
def thickness(C):
    U, S, Vh = torch.linalg.svd(C); M = U[:, :2]; P = (M @ M.conj().T).reshape(3, 3, 3, 3)
    return min(torch.linalg.eigvalsh(torch.einsum('abcb->ac', P))[0].item(), torch.linalg.eigvalsh(torch.einsum('abad->bd', P))[0].item())
for k in range(STARTS):
    A = torch.randn(D, 2, dtype=torch.complex128); C = A @ A.conj().T; C = C / torch.linalg.norm(C)
    def f(xn):
        x = torch.tensor(xn, requires_grad=True); r = R(squeeze(C, x)); r.backward(); return r.item(), x.grad.numpy().copy()
    res = minimize(f, np.zeros(18), jac=True, method='L-BFGS-B', options={'maxiter': 3000, 'gtol': 1e-12})
    X = squeeze(C, torch.tensor(res.x))
    print(f"  start {k}: R -> {res.fun:+.6f} (1 - 2c = {1 - 2 * c:+.6f}); |TrX|^2/||X||^2 = "
          f"{(abs(torch.trace(X)) ** 2 / torch.linalg.norm(X) ** 2).item():.6f}; thickness {thickness(X):.3f}; max|H| {np.abs(res.x).max():.2f}")
