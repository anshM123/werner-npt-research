# Moment (dual) test of SOS for f(x,y) = q_3(xx^T + yy^T) restricted to R^{d1}(x)R^{d2}(x)R^{d3}, with SCS.
# min L(f) s.t. moment matrix M_L (deg-2 monomials) PSD, L(|z|^4) = 1.  A clearly negative minimum => f is not SOS.
import sys, itertools, numpy as np, scipy.sparse as sps, cvxpy as cp
dims=tuple(int(a) for a in sys.argv[1].split(',')); D=int(np.prod(dims)); nv=2*D
mons=list(itertools.combinations_with_replacement(range(nv),2)); M=len(mons); midx={m:i for i,m in enumerate(mons)}
V=np.zeros((D,D,M))
for p in range(D):
    for q in range(D):
        V[p,q,midx[tuple(sorted((p,q)))]]+=1; V[p,q,midx[tuple(sorted((D+p,D+q)))]]+=1
n=len(dims); G0=np.zeros((M,M))
for r in range(n+1):
    for J in itertools.combinations(range(n),r):
        T=V.reshape(dims+dims+(M,))
        for s in sorted(J,reverse=True): T=np.trace(T,axis1=s,axis2=s+(T.ndim-1)//2)
        Tm=T.reshape(-1,M); G0+=(-0.5)**len(J)*(Tm.T@Tm)
quart={}; rows=[];cols=[]
for i,(a,b) in enumerate(mons):
    for j,(c,d) in enumerate(mons):
        k=quart.setdefault(tuple(sorted((a,b,c,d))),len(quart)); rows.append(i*M+j); cols.append(k)
B=sps.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(M*M,len(quart)))   # vec(M_L) = B y
fco=B.T@G0.reshape(-1)                                                          # L(f) = sum_alpha f_alpha y_alpha
Nn=np.zeros((M,M))
for a in range(nv):
    for b in range(nv): Nn[midx[(a,a)],midx[(b,b)]]+=1                          # |z|^4 = sum_ab z_a^2 z_b^2
nco=B.T@Nn.reshape(-1)
print(f'dims {dims}: {nv} vars, moment matrix {M}x{M}, {len(quart)} moments',flush=True)
y=cp.Variable(len(quart)); Mm=cp.reshape(B@y,(M,M),order='C')
pr=cp.Problem(cp.Minimize(fco@y),[(Mm+Mm.T)/2>>0, nco@y==1])
pr.solve(solver=cp.SCS,eps=1e-7,max_iters=200000,verbose=False)
print("status",pr.status," min L(f) =",pr.value,flush=True)
Mv=(B@y.value).reshape(M,M); Mv=(Mv+Mv.T)/2; ev=np.linalg.eigvalsh(Mv)
print("check: L(f) =",fco@y.value," L(|z|^4) =",nco@y.value," lambda_min(M_L) =",ev[0]," lambda_max =",ev[-1],flush=True)
