# Real factor Hessian of F_n(A,B) = 2^n q_n(AB^*) at the signed GHZ zero; usage: python ghz_zero_hessian.py nmax
# Real factor Hessian of F_n(A,B) = 2^n q_n(A B^*) at the signed GHZ zero C_n^0 = |0^n><0^n| - (-1)^n |1^n><1^n|.
# Claims (audit 71): PSD, kernel dimension 4*2^n + 4n + 2, smallest positive eigenvalue >= 4.
import torch, itertools, sys
torch.set_default_dtype(torch.float64)
def qn(C,n):
    tot=0
    for r in range(n+1):
        for J in itertools.combinations(range(n),r):
            T=C.reshape((3,)*(2*n))
            for s in sorted(J,reverse=True): T=torch.diagonal(T,dim1=s,dim2=s+T.dim()//2).sum(-1)
            m=int(round(T.numel()**0.5)); tot=tot+(-0.5)**len(J)*((T.reshape(m,m).real**2+T.reshape(m,m).imag**2).sum())
    return tot
for n in range(1,int(sys.argv[1])+1):
    D=3**n; sig=(-1)**n
    a=torch.zeros(D); a[0]=1; b=torch.zeros(D); b[sum(3**k for k in range(n))]=1
    A0=torch.stack([a,b],1).to(torch.complex128); B0=(A0@torch.diag(torch.tensor([1.0,-sig]).to(torch.complex128)))
    def F(x):
        A=torch.complex(x[0:2*D].reshape(D,2),x[2*D:4*D].reshape(D,2))+A0
        B=torch.complex(x[4*D:6*D].reshape(D,2),x[6*D:8*D].reshape(D,2))+B0
        return 2**n*qn(A@B.conj().T,n)
    x0=torch.zeros(8*D)
    print(f'n={n}: F at zero = {F(x0).item():.2e}',flush=True)
    H=torch.autograd.functional.hessian(F,x0)
    ev=torch.linalg.eigvalsh((H+H.T)/2)
    tol=1e-8; nz=int((ev.abs()<tol).sum()); pos=ev[ev>tol]
    print(f'  min eigenvalue {ev.min().item():.3e}, kernel dim {nz} (claim {4*2**n+4*n+2}), smallest positive {pos.min().item():.4f}',flush=True)
