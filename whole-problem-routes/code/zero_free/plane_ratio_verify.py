# Re-find the n=4 Schur violation (diagonal complex t), then verify in plain numpy and look for a rank-2 zero
# inside the polydisk: 0 in the numerical range of G(t) (4x4 Gram of the plane's rank-2 operators).
import numpy as np, torch, itertools, math, sys
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); d = 3
subsets = list(itertools.product([0,1], repeat=n))
def ptr(Xm, J):
    # partial trace of operator Xm (d^n x d^n) over sites in J (numpy, explicit)
    T = Xm.reshape([d]*(2*n)); m = n
    for i in sorted([i for i in range(n) if J[i]], reverse=True):
        T = np.trace(T, axis1=i, axis2=i+m); m -= 1
    k = n - sum(J)
    return T.reshape(d**k, d**k)
def Fpoly(B1, B2, t):
    # sum_J t^{|J|} <Tr_J B1, Tr_J B2>
    s = 0
    for J in subsets:
        s += t**sum(J) * np.vdot(ptr(B1, J), ptr(B2, J))
    return s
# torch search for diagonal t (one complex number) and plane maximising |X/D|
def pair(x, y, t):
    tx = x.reshape([d]*n); ty = y.reshape([d]*n); tot = 0
    for K in subsets:
        Kl = [i for i in range(n) if K[i]]; Kc = [i for i in range(n) if not K[i]]
        A = torch.tensordot(tx, ty.conj(), dims=(Kc, Kc)).reshape(d**len(Kl), d**len(Kl))
        B = torch.tensordot(ty, tx.conj(), dims=(Kc, Kc)).reshape(d**len(Kl), d**len(Kl))
        tot = tot + t**(n - len(Kl)) * torch.trace(A @ B)
    return tot
best = None
for r in range(4):
    torch.manual_seed(900 + r)
    U = torch.randn(d**n, 2, requires_grad=True); a = torch.tensor([0.0, 1.0], requires_grad=True)
    opt = torch.optim.Adam([U, a], lr=0.02)
    for it in range(1500):
        Q, _ = torch.linalg.qr(U); z = torch.complex(Q[:,0], Q[:,1])/math.sqrt(2)
        t = 0.5*torch.sigmoid(a[0])*torch.exp(1j*a[1])
        # D = sum_J t^|J| ||Tr_J zz*||^2 = sum_K t^{n-|K|} Tr[rho_K^2];  X = sum_K t^{n-|K|} Tr[rho_K rho_K^T]
        D = pair(z, z, t)
        # X: <Tr_J zz*, Tr_J zbar z^T> = Tr[Tr_{J}(zz*) Tr_J(zbar z^T)] ; use tensordot with zbar
        zb = z.conj()
        tx = z.reshape([d]*n); tb = zb.reshape([d]*n); X = 0
        for K in subsets:
            Kl = [i for i in range(n) if K[i]]; Kc = [i for i in range(n) if not K[i]]
            A = torch.tensordot(tx, tx.conj(), dims=(Kc, Kc)).reshape(d**len(Kl), d**len(Kl))
            B = torch.tensordot(tb, tb.conj(), dims=(Kc, Kc)).reshape(d**len(Kl), d**len(Kl))
            X = X + t**(n - len(Kl)) * torch.trace(A @ B)
        ratio = X.abs()/D.abs()
        opt.zero_grad(); (-ratio).backward(); opt.step()
    print(f"search run {r}: |X/D| = {ratio.item():.6f}, |t| = {abs(t.item()):.6f}, arg t = {np.angle(t.item()):.4f}", flush=True)
    if best is None or ratio.item() > best[0]:
        best = (ratio.item(), z.detach().numpy().copy(), complex(t.item()))
ratio, z, t = best
# shrink t slightly inside the open disk and verify in numpy
t = t * 0.999
u = np.sqrt(2)*z.real; v = np.sqrt(2)*z.imag
print("orthonormality check", np.dot(u,u), np.dot(v,v), np.dot(u,v))
Cz = np.outer(z, z.conj()); Czb = np.outer(z.conj(), z)
Dn = Fpoly(Cz, Cz, t); Xn = Fpoly(Cz, Czb, t)
print(f"numpy: |t| = {abs(t):.6f}, D = {Dn:.6f}, X = {Xn:.6f}, |X/D| = {abs(Xn/Dn):.6f}")
# Gram G(t) over basis of rank-<=2 operators in the plane: P, J, H1, H2
P = np.outer(u,u)+np.outer(v,v); Jm = np.outer(u,v)-np.outer(v,u); H1 = np.outer(u,u)-np.outer(v,v); H2 = np.outer(u,v)+np.outer(v,u)
basis = [P, Jm, H1, H2]
G = np.array([[Fpoly(Bi, Bj, t) for Bj in basis] for Bi in basis])
# F_C(t) = b^dagger G b for C = sum b_i B_i ; 0 in numerical range?  max_theta lambda_min(Herm(e^{i th} G))
ths = np.linspace(0, 2*np.pi, 3601)
vals = [np.min(np.linalg.eigvalsh((np.exp(1j*th)*G + np.conj(np.exp(1j*th)*G).T)/2)) for th in ths]
print(f"max_theta lambda_min(Re e^(i th) G) = {max(vals):.6e}   (negative => some rank-2 C in the plane has F_C(t) = 0 with |t| < 1/2)")
np.savez('plane_ratio_violation_n%d.npz' % n, z=z, t=t)
