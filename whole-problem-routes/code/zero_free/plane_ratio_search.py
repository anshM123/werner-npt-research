# Is r(t) = X(t)/D(t) a Schur function on the polydisk |t_s| <= 1/2 ?
# z isotropic (z^T z = 0), D = <zz|(x)(I + t_s F_s)|zz>,  X = <zbar z|(x)(I + t_s F_s)|z zbar>.
# Maximise |X|/|D| over z (u,v real orthonormal, z=(u+iv)/sqrt2) and complex t with |t_s|<=1/2.
import torch, itertools, sys, math
torch.set_default_dtype(torch.float64)
n = int(sys.argv[1]); runs = int(sys.argv[2]); real_t = len(sys.argv) > 3 and sys.argv[3] == 'real'
d = 3
subsets = list(itertools.product([0,1], repeat=n))
def pair_terms(x, y):
    # returns dict K -> <x y|F_K|x y> = Tr[rho_x^K rho_y^K]  (K = set where F acts)
    tx = x.reshape([d]*n); ty = y.reshape([d]*n); out = {}
    for K in subsets:
        Kl = [i for i in range(n) if K[i]]; Kc = [i for i in range(n) if not K[i]]
        A = torch.tensordot(tx, tx.conj(), dims=(Kc, Kc)).reshape(d**len(Kl), d**len(Kl))
        B = torch.tensordot(ty, ty.conj(), dims=(Kc, Kc)).reshape(d**len(Kl), d**len(Kl))
        out[K] = torch.trace(A @ B)
    return out
def cross_terms(z):
    # <zbar z| F_K |z zbar> = Tr[ Tr_{K^c}(z zbar^*)... ] compute directly: = sum over index contraction
    # <a b|F_K|c e> with a=zbar,b=z,c=z,e=zbar :  = Tr[ (Tr_{Kc} |c><a|) (Tr_{Kc} |e><b|) ]  with appropriate order
    tz = z.reshape([d]*n); tzb = z.conj().reshape([d]*n); out = {}
    for K in subsets:
        Kl = [i for i in range(n) if K[i]]; Kc = [i for i in range(n) if not K[i]]
        # M1 = Tr_{Kc} |z><zbar| : indices (K ket of z, K bra of zbar)
        M1 = torch.tensordot(tz, tzb.conj(), dims=(Kc, Kc)).reshape(d**len(Kl), d**len(Kl))
        M2 = torch.tensordot(tzb, tz.conj(), dims=(Kc, Kc)).reshape(d**len(Kl), d**len(Kl))
        out[K] = torch.trace(M1 @ M2)
    return out
best = 0
for r in range(runs):
    torch.manual_seed(31*n + r)
    U = torch.randn(d**n, 2, requires_grad=True)
    rad = torch.randn(n, requires_grad=True); ang = (2*math.pi*torch.rand(n)).requires_grad_()
    opt = torch.optim.Adam([U, rad, ang], lr=0.02)
    for it in range(2000):
        Q, _ = torch.linalg.qr(U)
        z = torch.complex(Q[:,0], Q[:,1])/math.sqrt(2)
        if real_t:
            t = 0.5*torch.tanh(rad).to(torch.complex128)
        else:
            t = 0.5*torch.sigmoid(rad)*torch.exp(1j*ang)
        Dt = pair_terms(z, z); Xt = cross_terms(z)
        D = 0; X = 0
        for K in subsets:
            w = 1
            for s in range(n):
                if K[s]: w = w * t[s]
            D = D + w*Dt[K]; X = X + w*Xt[K]
        ratio = X.abs()/D.abs()
        opt.zero_grad(); (-ratio).backward(); opt.step()
    best = max(best, ratio.item())
    print(f"run {r}: |X/D| = {ratio.item():.6f}  t = {[complex(round(v.real,3), round(v.imag,3)) for v in t.tolist()]}", flush=True)
print("max |X/D|", best)
