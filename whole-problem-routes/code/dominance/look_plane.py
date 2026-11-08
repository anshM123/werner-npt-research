import torch, sys
sys.argv=['dominance_test.py','5','1','g1','1']
src=open('dominance_test.py').read().split('best=-1e9')[0]; exec(src)
torch.manual_seed(1)
xr=torch.randn(3**n,dtype=dt,requires_grad=True); yr=torch.randn(3**n,dtype=dt,requires_grad=True)
opt=torch.optim.Adam([xr,yr],lr=0.02)
for it in range(400):
    x=xr/xr.norm(); y=yr-(yr@x)*x; y=y/y.norm()
    g,ep,hw,hs=gaps(x,y); opt.zero_grad(); (-g[1]).backward(); opt.step()
x=x.detach(); y=y.detach()
print('hist wedge',[round(v,4) for v in hw.tolist()]); print('hist sym  ',[round(v,4) for v in hs.tolist()])
for s in range(n):   # one-copy marginal spectra of the plane's projector and of x, y
    P=torch.outer(x,x)+torch.outer(y,y)
    T=P.reshape((3,)*n*2); perm=[s]+[t for t in range(n) if t!=s]; T=T.permute(perm+[p+n for p in perm]).reshape(3,3**(n-1),3,3**(n-1))
    m=torch.einsum('aibi->ab',T); print('copy',s,'marginal eigenvalues of P:',[round(v,4) for v in torch.linalg.eigvalsh(m).tolist()])
