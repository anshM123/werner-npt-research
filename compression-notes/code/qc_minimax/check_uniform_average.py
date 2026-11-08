# Check E_w Delta_s(w) = (3/4) q_n(C) - S/8 at t = 1/2 (uniform unit w in C^3), n = 2 and 3, random rank-2 C.
import numpy as np, itertools
rng=np.random.default_rng(1)
def qform(C,n):
    tot=0
    for r in range(n+1):
        for J in itertools.combinations(range(n),r):
            T=C.reshape((3,)*(2*n))
            keep=[s for s in range(n) if s not in J]
            # partial trace over J
            idx=list(range(2*n)); 
            for s in sorted(J,reverse=True):
                T=np.trace(T,axis1=s,axis2=s+T.ndim//2)
            m=3**len(keep); tot+=(-0.5)**len(J)*np.sum(np.abs(T.reshape(m,m))**2)
    return tot
for n in (2,3):
    D=3**n
    Y=rng.normal(size=(D,2))+1j*rng.normal(size=(D,2)); Z=rng.normal(size=(D,2))+1j*rng.normal(size=(D,2)); C=Y@Z.conj().T
    s=0; qn=qform(C,n)
    T=C.reshape(3,3**(n-1),3,3**(n-1)); S=sum(qform(T[i,:,j,:],n-1) for i in range(3) for j in range(3))
    vals=[]
    for _ in range(4000):
        w=rng.normal(size=3)+1j*rng.normal(size=3); w/=np.linalg.norm(w); P=np.eye(3)-np.outer(w,w.conj())
        PP=np.kron(P,np.eye(3**(n-1))); vals.append(qn-qform(PP@C@PP,n))
    print(f'n={n}: MC mean Delta = {np.mean(vals):.4f} +- {np.std(vals)/np.sqrt(len(vals)):.4f};  (3/4)q_n - S/8 = {0.75*qn-S/8:.4f}')
