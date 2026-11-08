"""Dense check of two operators that defeat canonical direction rules for Conjecture QC (NOTES.md, Section 1).
For each copy s, rho_s is the copy-s marginal of CC* + C*C; Delta_s(w) = q_n(C) - q_n(C_{s,w}), C_{s,w} = (P_s) C (P_s),
P = I - ww*.  Prints q_n/||C||^2 and Delta/||C||^2 at the three eigenvectors of rho_s at every copy.
usage: python check_rule_kills.py"""
import numpy as np, itertools
def q(C,n):
    tot=0
    for r in range(n+1):
        for J in itertools.combinations(range(n),r):
            T=C.reshape((3,)*(2*n))
            for s in sorted(J,reverse=True): T=np.trace(T,axis1=s,axis2=s+T.ndim//2)
            m=int(round(np.sqrt(T.size))); tot+=(-0.5)**len(J)*np.sum(np.abs(T.reshape(m,m))**2)
    return tot
def marg(X,n,s):
    T=X.reshape((3,)*(2*n))
    for t in sorted([u for u in range(n) if u!=s],reverse=True): T=np.trace(T,axis1=t,axis2=t+T.ndim//2)
    return T.reshape(3,3)
def compress(C,n,s,w):
    ops=[np.eye(3)]*n; ops[s]=np.eye(3)-np.outer(w,w.conj()); K=ops[0]
    for o in ops[1:]: K=np.kron(K,o)
    return K@C@K
for f,n in [('kill_n2_least_eigenvector.npy',2),('kill_n3_best_of_three_thinnest_copy.npy',3)]:
    C=np.load(f); nrm=np.sum(np.abs(C)**2); qn=q(C,n)
    print(f'{f}: n = {n}, rank {np.linalg.matrix_rank(C,tol=1e-9)}, q_n/|C|^2 = {qn/nrm:.5f}')
    for s in range(n):
        ev,V=np.linalg.eigh(marg(C@C.conj().T+C.conj().T@C,n,s))
        vals=[(qn-q(compress(C,n,s,V[:,k]),n))/nrm for k in range(3)]
        print(f'  copy {s}: relative least eigenvalue {ev[0]/ev.sum():.4f}; Delta at eigenvectors (least, middle, largest) = {np.round(vals,5)}')
print('(1 - sqrt 5)/8 =', (1-np.sqrt(5))/8)
