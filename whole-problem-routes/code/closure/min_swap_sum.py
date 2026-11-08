# Least value of <F_3> + <F_N> over filtered Werner forms G = (M (x) M')^* (I - F/2) (M (x) M') on C^3 (x) C^N, N = 3.
import torch, sys
torch.manual_seed(int(sys.argv[1]) if len(sys.argv)>1 else 0)
d=3; N=3; D=d*N
dt=torch.complex128
I=torch.eye(D*D,dtype=dt)
def swap(dims, which):
    # swap of the factors in 'which' between two replicas; space = (C^d x C^N)_1 (x) (C^d x C^N)_2
    n=len(dims); S=torch.zeros(D*D,D*D,dtype=dt)
    import itertools
    for idx in itertools.product(*[range(k) for k in dims+dims]):
        a=list(idx[:n]); b=list(idx[n:])
        a2=a[:]; b2=b[:]
        for w in which: a2[w],b2[w]=b[w],a[w]
        def lin(v):
            r=0
            for k,x in zip(dims+dims,v): r=r*k+x
            return r
        S[lin(a2+b2),lin(a+b)]=1
    return S
F3=swap([d,N],[0]); FN=swap([d,N],[1]); F=F3@FN
W=I-F/2
best=1e9
for run in range(6):
    Mr=torch.randn(D,D,dtype=torch.float64,requires_grad=True); Mi=torch.randn(D,D,dtype=torch.float64,requires_grad=True)
    Pr=torch.randn(D,D,dtype=torch.float64,requires_grad=True); Pi=torch.randn(D,D,dtype=torch.float64,requires_grad=True)
    opt=torch.optim.Adam([Mr,Mi,Pr,Pi],lr=0.03)
    for it in range(1500):
        M=torch.complex(Mr,Mi); Mp=torch.complex(Pr,Pi)
        K=torch.kron(M,Mp)
        G=K.conj().T@W@K
        val=(torch.trace(G@(F3+FN)).real)/torch.trace(G).real
        opt.zero_grad(); val.backward(); opt.step()
    v=val.item(); best=min(best,v)
    print('run',run,'min <F3>+<FN> =',round(v,6),flush=True)
print('best',best,' (target W3xW3 value = -0.4)')
