# Check E_U E_{s,s'} [(-3/5)^{Hamming(s,s')}] = (2/5)^n q_n(P)/4 for rank-two projectors P (s,s' iid from <s|UPU*|s>/2).
import numpy as np, itertools
from scipy.stats import unitary_group
rng=np.random.default_rng(3)
def q(C,n):
    tot=0
    for r in range(n+1):
        for J in itertools.combinations(range(n),r):
            T=C.reshape((3,)*(2*n))
            for s in sorted(J,reverse=True): T=np.trace(T,axis1=s,axis2=s+T.ndim//2)
            m=int(round(np.sqrt(T.size))); tot+=(-0.5)**len(J)*np.sum(np.abs(T.reshape(m,m))**2)
    return tot.real
for n in (2,3):
    D=3**n; X=rng.normal(size=(D,2))+1j*rng.normal(size=(D,2)); Qm,_=np.linalg.qr(X); P=Qm@Qm.conj().T
    K1=np.full((3,3),-0.6)+np.eye(3)*1.6; K=K1
    for _ in range(n-1): K=np.kron(K,K1)
    vals=[]
    for _ in range(20000):
        U=unitary_group.rvs(3,random_state=rng)
        for _ in range(n-1): U=np.kron(U,unitary_group.rvs(3,random_state=rng))
        p=np.real(np.einsum('ij,jk,ik->i',U,P,U.conj()))/2
        vals.append(p@K@p)
    print(f'n={n}: MC {np.mean(vals):.5f} +- {np.std(vals)/np.sqrt(len(vals)):.5f}   formula (2/5)^n q_n(P)/4 = {(0.4**n)*q(P,n)/4:.5f}')
