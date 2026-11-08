# Stochastic dominance test for the swap-test form of the endpoint.  Real plane V = span{x,y} (orthonormal) on n qutrit
# copies.  R = set of copies where the local swap measurement {P_S,P_A} on V (x) V gives 'antisymmetric'.
#   mu_wedge(R): distribution for x^y,   mu_sym(R): average over Sym^2 V (maximally mixed, 3-dim).
# EP_n  <=>  E_wedge[3^|R|] <= 3 E_sym[3^|R|].   Dominance D_n: P_wedge(|R| >= k+1) <= P_sym(|R| >= k) for all k.
import torch, itertools, sys
n=int(sys.argv[1]); runs=int(sys.argv[2]); which=sys.argv[3]   # which in {g1,g2,ep}
torch.manual_seed(int(sys.argv[4]) if len(sys.argv)>4 else 0)
dt=torch.float64; subs=[K for r in range(n+1) for K in itertools.combinations(range(n),r)]
def mats(v,K):
    Kc=[s for s in range(n) if s not in K]
    t=v.reshape((3,)*n).permute(list(K)+Kc) if n>0 else v
    return t.reshape(3**len(K),3**len(Kc))
def dists(x,y):
    f={'xx':[],'yy':[],'w':[],'v':[]}
    for K in subs:
        Mx,My=mats(x,K),mats(y,K)
        pxy=((Mx.T@My)**2).sum()            # Tr[rho^x_K rho^y_K]
        cr=((Mx@My.T)**2).sum()             # ||Tr_{K^c}(x y^T)||^2
        f['xx'].append(((Mx.T@Mx)**2).sum()); f['yy'].append(((My.T@My)**2).sum())
        f['w'].append(pxy-cr); f['v'].append(pxy+cr)
    P={}
    for key in f:
        fk=torch.stack(f[key]); pr=[]
        for R in subs:
            sg=torch.tensor([(-1.0)**len(set(K)&set(R)) for K in subs],dtype=dt)
            pr.append((sg*fk).sum()/2**n)
        P[key]=torch.stack(pr)
    size=torch.tensor([len(R) for R in subs])
    hist=lambda p: torch.stack([p[size==k].sum() for k in range(n+1)])
    hw=hist(P['w']); hs=(hist(P['xx'])+hist(P['yy'])+hist(P['v']))/3
    return hw,hs
def gaps(x,y):
    hw,hs=dists(x,y)
    tailw=lambda k: hw[k:].sum(); tails=lambda k: hs[k:].sum()
    g=[tailw(k+1)-tails(k) for k in range(n)]
    ep=(hw*3.0**torch.arange(n+1)).sum()-3*(hs*3.0**torch.arange(n+1)).sum()
    return torch.stack(g),ep,hw,hs
best=-1e9
for run in range(runs):
    xr=torch.randn(3**n,dtype=dt,requires_grad=True); yr=torch.randn(3**n,dtype=dt,requires_grad=True)
    opt=torch.optim.Adam([xr,yr],lr=0.02)
    for it in range(400):
        x=xr/xr.norm(); y=yr-(yr@x)*x; y=y/y.norm()
        g,ep,hw,hs=gaps(x,y)
        obj = ep if which=='ep' else (g[1] if which=='g1' else g[3])
        opt.zero_grad(); (-obj).backward(); opt.step()
    best=max(best,obj.item())
    print(f'run {run}: {which} = {obj.item():+.6f} | all gaps {[round(v,5) for v in g.tolist()]} | EP value {ep.item():+.5f}',flush=True)
print('max',which,best)
