// The energy model against the relaxed form of ../typeI_G_m1_relax (Lemma R1), entry by entry (NOTES.md, Section 4).
// Run it with `python run_js.py relax_comparison.js sfQReviewA510RelaxAllVerify`.

function sfGcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b){const t=a%b;a=b;b=t;}return a;}
class SFQ{
 constructor(a,b=1n){if(a instanceof SFQ){this.n=a.n;this.d=a.d;return;}a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero denominator");if(b<0n){a=-a;b=-b;}const g=sfGcd(a,b);this.n=a/g;this.d=b/g;}
 add(y){y=sfq(y);return new SFQ(this.n*y.d+y.n*this.d,this.d*y.d);}
 sub(y){y=sfq(y);return new SFQ(this.n*y.d-y.n*this.d,this.d*y.d);}
 mul(y){y=sfq(y);return new SFQ(this.n*y.n,this.d*y.d);}
 div(y){y=sfq(y);return new SFQ(this.n*y.d,this.d*y.n);}
 neg(){return new SFQ(-this.n,this.d);}
 sq(){return this.mul(this);}
 eq(y){y=sfq(y);return this.n===y.n&&this.d===y.d;}
 gt(y){y=sfq(y);return this.n*y.d>y.n*this.d;}
 ge(y){y=sfq(y);return this.n*y.d>=y.n*this.d;}
 str(){return this.d===1n?String(this.n):this.n+"/"+this.d;}
 num(){const a=Number(this.n),b=Number(this.d);if(Number.isFinite(a)&&Number.isFinite(b))return a/b;const s=this.n<0n?-1:1,an=(this.n<0n?-this.n:this.n).toString(),ad=this.d.toString(),k=16;return s*(Number(an.slice(0,k))/10**(Math.min(k,an.length)-1))/(Number(ad.slice(0,k))/10**(Math.min(k,ad.length)-1))*10**(an.length-ad.length);}
}
const sfq=(a,b)=>a instanceof SFQ?a:new SFQ(a,b);
const sfsum=a=>a.reduce((t,x)=>t.add(x),sfq(0));
const sfdot=(a,b)=>sfsum(a.map((x,i)=>x.mul(b[i])));
function sfdet(M){const a=M.map(r=>r.map(x=>sfq(x)));let v=sfq(1);for(let p=0;p<a.length;p++){let r=p;while(r<a.length&&a[r][p].eq(0))r++;if(r===a.length)return sfq(0);if(r!==p){[a[r],a[p]]=[a[p],a[r]];v=v.neg();}let d=a[p][p];v=v.mul(d);for(let i=p+1;i<a.length;i++){const f=a[i][p].div(d);for(let j=p+1;j<a.length;j++)a[i][j]=a[i][j].sub(f.mul(a[p][j]));}}return v;}
function sfInvE(w,K){const d=w.map(x=>K.add(x.sq())),eta=w.map((x,i)=>x.div(d[i])),G=sfq(2).add(sfsum(w.map((x,i)=>x.sq().div(d[i]))));return w.map((_,i)=>w.map((_,j)=>(i===j?sfq(1).div(d[i]):sfq(0)).sub(eta[i].mul(eta[j]).div(G))));}
function sfMatVec(A,v){return A.map(r=>sfdot(r,v));}
function sfT1(w,r,v){const z=v.map(x=>x.mul(w[r]).div(2));z[r]=z[r].add(sfdot(w,v).div(2)).add(w[r].mul(v[r]));return z;}
function sfCUForm(w,sigma,X){const S=sfsum(w.map(x=>x.sq())),K=sigma.add(S.div(2)),Ei=sfInvE(w,K);const tr=sfsum(X.map((r,i)=>r[i])),beta=sfdot(w,sfMatVec(X,w));let f=sfsum(X.flat().map(x=>x.sq())).div(2).add(sfsum(X.map((r,i)=>r[i].sq()))).add(tr.sq().div(2)).sub(tr.sub(beta.div(K.mul(2))).sq().div(2));const V=X.map(v=>w.map((_,r)=>sfT1(w,r,v)));for(let r=0;r<w.length;r++)for(let h=0;h<w.length;h++)f=f.sub(sfdot(V[r][h],sfMatVec(Ei,V[h][r])));return f;}
function sfHardForm(w,sigma,X){const S=sfsum(w.map(x=>x.sq())),K=sigma.add(S.div(2)),rho=w.map(x=>x.sq().div(K.add(x.sq()))),G=sfq(2).add(sfsum(rho));let f=sfq(0);for(let r=0;r<w.length;r++)for(let h=0;h<w.length;h++){if(r===h)f=f.add(sfq(3,2).sub(rho[r]).mul(X[r][r].sq()));else f=f.add(X[r][h].sq().div(2));f=f.add(rho[r].mul(rho[h]).mul(X[r][h]).mul(X[h][r]).div(G));}return f;}
function sfHardProject(w,Z){const S=sfdot(w,w),Pi=w.map((a,i)=>w.map((b,j)=>(i===j?sfq(1):sfq(0)).sub(a.mul(b).div(S))));const m=(A,B)=>A.map(r=>B[0].map((_,j)=>sfsum(r.map((x,i)=>x.mul(B[i][j])))));return m(m(Pi,Z),Pi);}

const sfQT=A=>A[0].map((_,j)=>A.map(r=>r[j]));
const sfQM=(A,B)=>A.map(r=>B[0].map((_,j)=>sfsum(r.map((x,k)=>x.mul(B[k][j])))));
function sfQInv(A){const n=A.length,D=A.map((r,i)=>[...r,...r.map((_,j)=>sfq(i===j?1:0))]);for(let p=0;p<n;p++){const v=D[p][p];if(v.eq(0))throw Error("zero inverse pivot");for(let j=p;j<2*n;j++)D[p][j]=D[p][j].div(v);for(let i=0;i<n;i++)if(i!==p){let z=D[i][p];for(let j=p;j<2*n;j++)D[i][j]=D[i][j].sub(z.mul(D[p][j]));}}return D.map(r=>r.slice(n));}
function sfQLDL(A){const n=A.length,D=A.map(r=>r.slice()),piv=[];for(let p=0;p<n;p++){let v=D[p][p];piv.push(v);if(v.eq(0))return piv;for(let i=p+1;i<n;i++)for(let j=p+1;j<n;j++)D[i][j]=D[i][j].sub(D[i][p].mul(D[p][j]).div(v));}return piv;}

function sfQInertia(M){
 let A=M.map(row=>row.slice()),positive=0,negative=0;
 while(A.length){
  const n=A.length;let i=A.findIndex((row,j)=>!row[j].eq(0));
  if(i>=0){if(i!==0){[A[0],A[i]]=[A[i],A[0]];for(const row of A)[row[0],row[i]]=[row[i],row[0]];}
   const p=A[0][0];if(p.gt(0))positive++;else negative++;
   A=A.slice(1).map((row,j)=>row.slice(1).map((v,k)=>v.sub(row[0].mul(A[0][k+1]).div(p))));continue;
  }
  let ij=null;for(let r=0;r<n&&!ij;r++)for(let c=r+1;c<n;c++)if(!A[r][c].eq(0)){ij=[r,c];break;}
  if(!ij)return{positive,negative,zero:n,rank:positive+negative};
  const [r,c]=ij,order=[r,c,...Array.from({length:n},(_,i)=>i).filter(i=>i!==r&&i!==c)];A=order.map(i=>order.map(j=>A[i][j]));const b=A[0][1];positive++;negative++;
  A=A.slice(2).map((row,j)=>row.slice(2).map((v,k)=>v.sub(row[0].mul(A[1][k+2]).add(row[1].mul(A[0][k+2])).div(b))));
 }
 return{positive,negative,zero:0,rank:positive+negative};
}

const spAdd=(a,b)=>{const c=Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]??sfq(0)).add(b[i]??sfq(0)));while(c.length>1&&c.at(-1).eq(0))c.pop();return c;};
const spNeg=a=>a.map(v=>v.neg()),spSub=(a,b)=>spAdd(a,spNeg(b)),spC=v=>[sfq(v)],spScale=(a,v)=>a.map(x=>x.mul(v));
function spMul(a,b){const c=Array(a.length+b.length-1).fill(sfq(0));for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]=c[i+j].add(a[i].mul(b[j]));while(c.length>1&&c.at(-1).eq(0))c.pop();return c;}
const spPow=(a,n)=>{let z=spC(1);while(n){if(n&1)z=spMul(z,a);n>>=1;if(n)a=spMul(a,a);}return z;};
const spEval=(a,x)=>{let z=sfq(0);for(let i=a.length-1;i>=0;i--)z=z.mul(x).add(a[i]);return z;};
const spProd=(...args)=>args.reduce(spMul,spC(1)),spSum=args=>args.reduce(spAdd,spC(0));

function sfQMeasureModelGeneralExpanded(wt,x,mu,A,C,Z,rho,coupling=sfq(5,2)){
 const m=wt.length,sigma=sfdot(wt,wt),K=sigma.add(sfq(1,2)),avg=v=>sfdot(mu,v),p=avg(C),F=avg(x.map(s=>sfq(1).div(K.add(s)))),G=sfq(2).add(F),Af=avg(A.map((v,i)=>v.div(K.add(x[i])))),Cf=avg(C.map((v,i)=>v.div(K.add(x[i])))),IA=avg(A.map((v,i)=>v.sq().div(K.add(x[i])))).sub(Af.sq().div(G)),IC=avg(C.map((v,i)=>v.sq().div(K.add(x[i])))).sub(Cf.sq().div(G));
 let base=sigma.mul(avg(A.map((v,i)=>v.sq().add(C[i].sq()))).sub(sigma.mul(IA.add(IC))).sub(sigma.mul(p.sq()).div(K.sq().mul(2))).add(p.mul(Af).mul(2).div(G)).add(avg(x.map((s,i)=>s.mul(A[i]).mul(C[i]).div(K.add(s)))).mul(4).div(G))).add(sfdot(rho,rho).mul(sfq(3,2))).add(sfsum(rho).sq().div(2)).add(p.mul(sfsum(rho))).add(sfsum(Z.map(z=>avg(z.map(v=>v.sq())))));
 function invB(W,k,vm,vf,um,uf){const ds=W.map(w=>k.add(w.sq())),d=x.map(s=>k.add(s)),L=sfq(2).add(sfsum(W.map((w,j)=>w.sq().div(ds[j])))).add(avg(d.map(v=>sfq(1).div(v)))),left=sfsum(W.map((w,j)=>w.mul(vm[j]).div(ds[j]))).add(avg(vf.map((v,i)=>v.div(d[i])))),right=sfsum(W.map((w,j)=>w.mul(um[j]).div(ds[j]))).add(avg(uf.map((v,i)=>v.div(d[i]))));return sfsum(vm.map((v,j)=>v.mul(um[j]).div(ds[j]))).add(avg(vf.map((v,i)=>v.mul(uf[i]).div(d[i])))).sub(left.mul(right).div(L));}
 const ell=[],betas=[],poles=[],denominators=[G];let Ly=sfq(0);
 for(let j=0;j<m;j++){
 const st=wt[j].sq(),k=K.sub(st.div(2)),z=avg(Z[j]),beta=wt[j].mul(rho[j]).add(z),bf=x.map((s,i)=>beta.div(2).add(s.mul(Z[j][i]))),W=[wt[j]],chiA=invB(W,k,[sfq(1,2)],x.map(()=>wt[j].div(-2)),[sfq(0)],A.map((a,i)=>a.mul(bf[i]))),chiC=invB(W,k,[p.div(2)],C.map(c=>c.mul(wt[j]).div(-2)),[rho[j].div(2)],bf.map((b,i)=>b.add(Z[j][i].div(2))));
 base=base.sub(chiA.add(chiC).mul(2));Ly=Ly.sub(K.mul(2).mul(rho[j]).add(wt[j].mul(z)).div(K.mul(2).sub(st)));
 const dm=k.add(st),d=x.map(s=>k.add(s)),L=sfq(2).add(st.div(dm)).add(avg(d.map(v=>sfq(1).div(v)))),ym=wt[j].mul(rho[j]).add(z.div(2)),yf=Z[j].map(v=>v.mul(wt[j]).div(-2)),moment=wt[j].mul(ym).div(dm).add(avg(yf.map((v,i)=>v.div(d[i])))),em=ym.sub(wt[j].mul(moment).div(L)).div(dm),ef=yf.map((v,i)=>v.sub(moment.div(L)).div(d[i]));
 ell.push(Z[j].map((v,i)=>v.neg().sub(em).add(wt[j].mul(ef[i]))));betas.push(beta);poles.push(k);denominators.push(L);
 }
 const LX=sigma.mul(p).neg().div(K);base=base.sub(LX.mul(Ly)).sub(Ly.sq().div(2)).sub(sfsum(ell.map(z=>avg(z.map(v=>v.sq())))).div(4)).add(avg(x.map((_,i)=>sfsum(wt.map((w,j)=>w.mul(ell[j][i]))).sq())).div(sigma.mul(4)));
 for(let j=0;j<m;j++)for(let h=j+1;h<m;h++){
 const k=K.sub(wt[j].sq().add(wt[h].sq()).div(2)),W=[wt[j],wt[h]],v1m=[wt[h].mul(rho[j]).div(-2),betas[j].div(2)],v2m=[betas[h].div(2),wt[j].mul(rho[h]).div(-2)],v1=Z[j].map(v=>v.mul(wt[h]).div(-2)),v2=Z[h].map(v=>v.mul(wt[j]).div(-2));
 base=base.sub(invB(W,k,v1m,v1,v2m,v2).mul(2));denominators.push(sfq(2).add(sfsum(W.map(w=>w.sq().div(k.add(w.sq()))))).add(avg(x.map(s=>sfq(1).div(k.add(s))))));
 }
 const V=x.map((s,i)=>sigma.mul(A[i].add(C[i]).div(K.add(s)).sub(p.div(K.mul(K.add(s)).mul(2))))),PP=[x.map(s=>sigma.div(G.mul(K.add(s)))),x.map((s,i)=>sigma.div(G).mul(s.mul(C[i]).div(K.add(s)).add(p.div(K.add(s).mul(2)))))],QQ=[x.map((s,i)=>s.mul(A[i]).div(K.add(s))),x.map(s=>sfq(1).div(K.add(s)))];
 for(let j=0;j<m;j++){const k=poles[j],L=denominators[j+1],b=wt[j].sq().mul(rho[j]).add(wt[j].mul(avg(Z[j])));for(let i=0;i<x.length;i++)V[i]=V[i].add(wt[j].mul(Z[j][i]).div(k.add(x[i]).mul(2))).sub(b.div(k.mul(k.add(x[i])).mul(4)));
 PP.push(x.map((s,i)=>wt[j].mul(s).mul(Z[j][i]).div(k.add(s).mul(2)).add(b.div(k.add(s).mul(4))).div(L)));QQ.push(x.map(s=>sfq(1).div(k.add(s))));
 }
 const center=v=>{const v0=avg(v);return v.map(v=>v.sub(v0));},P=PP.map(center),Q=QQ.map(center);
 let norm=avg(x.map((s,i)=>s.mul(sfq(1).sub(s.mul(2))).mul(V[i].sq()))).add(avg(x.map((s,i)=>s.mul(V[i]))).sq());
 for(let j=0;j<P.length;j++)norm=norm.sub(avg(x.map((s,i)=>s.mul(V[i]).mul(P[j][i]).mul(Q[j][i]))).mul(2));
 for(let j=0;j<P.length;j++)for(let h=0;h<P.length;h++)norm=norm.add(avg(P[j].map((v,i)=>v.mul(P[h][i]))).mul(avg(Q[j].map((v,i)=>v.mul(Q[h][i])))));
 return{base,norm,model:base.sub(coupling.mul(norm)),denominators,poles,PP,QQ,P,Q,V};
}

function sfQTriangleRemainder(wt,x,mu,A,C,Z,rho){
 const r=sfQMeasureModelGeneralExpanded(wt,x,mu,A,C,Z,rho),avg=v=>sfdot(mu,v),B=avg(r.V.map((v,i)=>x[i].add(sfq(1,4)).mul(v.sq()))).add(avg(r.V).sq().div(4));
 let R=sfq(0);for(let i=0;i<r.P.length;i++)for(let j=0;j<r.P.length;j++)R=R.add(avg(r.P[i].map((v,k)=>v.mul(r.P[j][k]))).mul(avg(r.Q[i].map((v,k)=>v.mul(r.Q[j][k])))));
 return{base:r.base,B,R,remainder:r.base.sub(B).sub(R.mul(8))};
}
function sfQTriangleRemainderMatrix(wt,x,mu){
 const m=wt.length,u=x.length,n=(m+2)*u+m-1,form=v=>{
 let A=Array(u).fill(sfq(0));for(let j=0;j<u-1;j++){A[j]=v[j];A[u-1]=A[u-1].sub(mu[j].mul(v[j]).div(mu[u-1]));}
 const C=v.slice(u-1,2*u-1),Z=Array.from({length:m},(_,j)=>v.slice(2*u-1+j*u,2*u-1+(j+1)*u)),rho=v.slice((m+2)*u-1);return sfQTriangleRemainder(wt,x,mu,A,C,Z,rho).remainder;
 };
 const basis=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>sfq(i===j?1:0))),val=basis.map(form),M=val.map((v,i)=>val.map((_,j)=>i===j?v:sfq(0)));
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)M[i][j]=M[j][i]=form(basis[i].map((v,k)=>v.add(basis[j][k]))).sub(val[i]).sub(val[j]).div(2);
 return M;
}


const sfQZero=n=>Array.from({length:n},()=>Array.from({length:n},()=>sfq(0)));
const sfQDiag=v=>v.map((x,i)=>v.map((_,j)=>i===j?x:sfq(0)));
const sfQSub=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.sub(B[i][j])));
const sfQAdd=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.add(B[i][j])));
const sfQScale=(A,c)=>A.map(row=>row.map(x=>x.mul(c)));
const sfQCong=(G,A)=>sfQM(sfQT(G),sfQM(A,G));
function sfQP6FormMatrix(n,fn){const b=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>sfq(i===j?1:0))),v=b.map(fn),M=sfQDiag(v);for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)M[i][j]=M[j][i]=fn(b[i].map((z,k)=>z.add(b[j][k]))).sub(v[i]).sub(v[j]).div(2);return M;}

function sfQPCommonData(wt,x,mu,A,C0,Z,rho){
 const t=sfdot(wt,wt),K=t.add(sfq(1,2)),f=x.map(s=>sfq(1).div(K.add(s))),F=sfdot(mu,f),alpha=sfq(1).div(sfq(2).add(F)),cp=t.sub(t.sq().mul(alpha).mul(F).mul(2)).sub(t.sq().div(K.sq().mul(2))),b0=t.sq().mul(alpha).mul(sfdot(mu,A.map((z,i)=>z.add(C0[i]).mul(f[i])))).mul(-2),b=[];
 for(let j=0;j<wt.length;j++){
 const a=wt[j].sq(),w=wt[j],k=K.sub(a.div(2)),delta=k.add(a),D0=sfq(2).add(a.div(delta)),Cstar=sfq(2).add(a.add(1).div(delta)),v=k.sub(sfq(1,2)),g=x.map(s=>sfq(1).div(k.add(s))),B=sfdot(mu,g),gam=sfq(1).div(D0.add(B)),h=sfdot(mu,Z[j]),Zg=sfdot(mu,Z[j].map((z,i)=>z.mul(g[i]))),bR=sfq(1).sub(sfq(1).div(delta.mul(2))).sub(t.div(k)).add(gam.mul(a).div(2).mul(D0.mul(B).add(sfq(1).div(delta.sq())))).div(2),bH=w.div(4).mul(Cstar.mul(sfq(1).sub(gam.mul(a).div(delta))).sub(t.div(K.mul(k)))),bg=w.mul(gam).mul(Cstar).mul(v).div(-2);
 b.push(bR.mul(rho[j]).add(bH.mul(h)).add(bg.mul(Zg)));
 }
 const lb=t.mul(t.sq().add(t).add(sfq(1,2))).div(K.sq().mul(t.add(1)).mul(2));if(!cp.ge(lb)||!lb.gt(0))throw Error("common p coefficient lower bound");
 return{t,K,F,alpha,cp,b0,b,source:b0.add(sfsum(b)),lower:lb};
}
function sfQPCommonIdentityGuard(wt,x,mu){
 const m=wt.length,u=x.length,n=2*(u-1)+m*u+m,field=v=>[...v,sfdot(mu.slice(0,u-1),v).neg().div(mu[u-1])],form=z=>{
 const A=field(z.slice(0,u-1)),C0=field(z.slice(u-1,2*(u-1))),Z=Array.from({length:m},(_,j)=>z.slice(2*(u-1)+j*u,2*(u-1)+(j+1)*u)),rho=z.slice(2*(u-1)+m*u),d=sfQPCommonData(wt,x,mu,A,C0,Z,rho),at=p=>sfQMeasureModelGeneralExpanded(wt,x,mu,A,C0.map(v=>v.add(p)),Z,rho).base,b=at(sfq(1)).sub(at(sfq(-1))).div(4),c=at(sfq(1)).add(at(sfq(-1))).div(2).sub(at(sfq(0)));
 if(!b.eq(d.source)||!c.eq(d.cp))throw Error("all-m common p identity");return true;
 },basis=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>sfq(i===j?1:0)));for(const v of basis)form(v);return{m,u,linearCoefficientChecks:n,pivotChecks:n};
}
function sfQPurePSource(wU,t){
 const x=wU.map(w=>w.sq()),mu=x,u=x.length,K=t.add(sfq(1,2)),f=x.map(s=>sfq(1).div(K.add(s))),F=sfdot(mu,f),alpha=sfq(1).div(sfq(2).add(F)),Pi=wU.map((w,i)=>wU.map((v,j)=>sfq(i===j?1:0).sub(w.mul(v)))),v=wU.map((w,i)=>w.mul(f[i].sub(F))),diag=wU.map((_,i)=>wU.map((_,j)=>i===j?x[i].mul(f[i]):sfq(0))),PD=sfQM(sfQM(Pi,diag),Pi),R=PD.map((row,i)=>row.map((z,j)=>t.sq().mul(z.div(K).add(alpha.mul(v[i]).mul(v[j])))));
 if(!sfdot(wU,wU).eq(1))throw Error("normalize free weights");
 return{R,x,F,alpha,K,PD};
}
function sfQPurePSourceGuard(wU,wT){
 const t=wT.sq(),u=wU.length,{R,x,F,alpha,K,PD}=sfQPurePSource(wU,t),zero=x.map(()=>sfq(0)),direct=sfQTriangleRemainder([wT],x,x,zero,x.map(()=>sfq(1)),[zero],[sfq(0)]),norm=sfsum(R.flat().map(z=>z.sq())),sourceN=sfQMeasureModelGeneralExpanded([wT],x,x,zero,x.map(()=>sfq(1)),[zero],[sfq(0)]).norm,cp=t.sub(t.sq().mul(alpha).mul(F).mul(2)).sub(t.sq().div(K.sq().mul(2))),pivots=sfQInertia(R),gapHalf=R.map((row,i)=>row.map((z,j)=>sfq(i===j?1:0).div(2).sub(z)));
 const compression=sfQInertia(PD.map((row,i)=>row.map((z,j)=>sfq(i===j?1:0).div(K.mul(2).add(1)).sub(z))));if(compression.negative)throw Error("compressed resolvent multiplier");
 if(!norm.eq(sourceN)||!direct.base.eq(cp)||!direct.remainder.gt(0)||sfQInertia(gapHalf).negative||pivots.negative)throw Error("pure p source guard");
 if(!cp.mul(sfq(29,72)).gt(norm))throw Error("pure p source relative bound");
 // Use a rational full basis of matrices with X w_U=w_U^T X=0, without requiring orthonormal square roots.
 const H=Array.from({length:u},(_,i)=>Array.from({length:u-1},(_,j)=>i===j?sfq(1):i===u-1?wU[j].neg().div(wU[u-1]):sfq(0))),n=(u-1)**2;
 const form=z=>{const V=Array.from({length:u-1},(_,i)=>z.slice(i*(u-1),(i+1)*(u-1))),X=sfQM(sfQM(H,V),sfQT(H)),rp=sfsum(X.flatMap((row,i)=>row.map((v,j)=>v.mul(R[i][j])))),fn=sfsum(X.flat().map(v=>v.sq()));return sfHardForm(wU,t,X).sub(rp.sq().div(cp)).sub(fn.div(72));},M=sfQP6FormMatrix(n,form),I=sfQInertia(M);
 if(I.negative||I.zero)throw Error("common p completed hard coercivity");
 return{sigma:t.str(),u,sourcePSD:pivots,compression,halfOperatorMargin:sfQInertia(gapHalf),relativeFrobeniusRatio:norm.div(cp).str(),pureEnergy:direct.remainder.str(),hardQuotientDimension:n,hardMargin:I};
}


function sfQSoftCommonSources(wt,x,mu,Z,rho){
 const zero=x.map(()=>sfq(0)),d=sfQPCommonData(wt,x,mu,zero,zero,Z,rho),{t,K,alpha}=d,JL=x.map(()=>sfq(0)),JD=x.map(()=>sfq(0));
 for(let j=0;j<wt.length;j++){
 const a=wt[j].sq(),w=wt[j],k=K.sub(a.div(2)),delta=k.add(a),D0=sfq(2).add(a.div(delta)),C=sfq(2).add(a.add(1).div(delta)),v=k.sub(sfq(1,2)),g=x.map(s=>sfq(1).div(k.add(s))),B=sfdot(mu,g),gam=sfq(1).div(D0.add(B)),h=sfdot(mu,Z[j]),ZZ=Z[j].map(z=>z.sub(h)),r=w.mul(rho[j]).sub(t.mul(2).sub(a).mul(h)),mZ=sfdot(mu,ZZ.map((z,i)=>z.mul(g[i]))),env=t.sub(a);
 for(let i=0;i<x.length;i++){const s=x[i],half=w.mul(g[i]).div(2),brL=gam.mul(t.mul(2).add(a.mul(2)).add(1)).div(delta).mul(r).sub(gam.mul(env).div(delta).mul(h)).add(gam.mul(v).mul(mZ)),brD=gam.div(delta).mul(r).add(gam.mul(env).div(delta).mul(h)).sub(gam.mul(v).mul(mZ));
 JL[i]=JL[i].add(half.mul(gam.mul(C).mul(s).add(s.add(sfq(1,2))).mul(ZZ[i]).add(brL)));
 JD[i]=JD[i].add(half.mul(gam.mul(C).mul(s).sub(s.add(sfq(1,2))).mul(ZZ[i]).add(brD)));}
 }
 return{...d,b:sfsum(d.b),JL,JD};
}
function sfQSoftCommonCovariance(m,f,J,mu,c){
 const W=sfdot(mu,m.map(z=>sfq(1).div(z))),fmean=sfdot(mu,f.map((z,i)=>z.div(m[i]))).div(W),Jmean=sfdot(mu,J.map((z,i)=>z.div(m[i]))).div(W),uf=f.map((z,i)=>z.sub(fmean).div(m[i])),uJ=J.map((z,i)=>z.sub(Jmean).div(m[i])),V=sfdot(mu,f.map((z,i)=>z.mul(uf[i]))),U=sfdot(mu,J.map((z,i)=>z.mul(uf[i]))),JJ=sfdot(mu,J.map((z,i)=>z.mul(uJ[i]))),den=sfq(1).add(c.mul(V));
 if(!den.gt(0))throw Error("soft inverse rank denominator");
 const optimizer=uJ.map((z,i)=>z.neg().add(c.mul(U).div(den).mul(uf[i]))),benefit=JJ.sub(c.mul(U.sq()).div(den));
 if(!sfdot(mu,optimizer).eq(0))throw Error("soft inverse mean constraint");return{W,V,U,den,benefit,optimizer};
}
function sfQSoftCommonInverse(wt,x,mu,Z,rho){
 const d=sfQSoftCommonSources(wt,x,mu,Z,rho),{t,K,alpha,cp,b,JL,JD}=d,f=x.map(s=>sfq(1).div(K.add(s))),mL=x.map((s,i)=>s.add(sfq(1,2)).add(alpha.mul(s).mul(2)).mul(f[i]).mul(t).mul(2)),mD=x.map((s,i)=>s.add(sfq(1,2)).sub(alpha.mul(s).mul(2)).mul(f[i]).mul(t).mul(2)),cD=t.sq().mul(alpha).mul(2),cL=cD.sub(t.sq().sq().mul(alpha.sq()).mul(16).div(cp)),shift=t.sq().mul(alpha).mul(b).mul(4).div(cp),L=sfQSoftCommonCovariance(mL,f,JL.map((z,i)=>z.add(shift.mul(f[i]))),mu,cL),D=sfQSoftCommonCovariance(mD,f,JD,mu,cD),p=t.sq().mul(alpha).mul(sfdot(mu,L.optimizer.map((z,i)=>z.mul(f[i])))).mul(4).sub(b).div(cp),A=L.optimizer.map((z,i)=>z.add(D.optimizer[i])),C=L.optimizer.map((z,i)=>z.sub(D.optimizer[i]).add(p)),zero=x.map(()=>sfq(0)),letter=sfQMeasureModelGeneralExpanded(wt,x,mu,zero,zero,Z,rho).base,value=letter.sub(b.sq().div(cp)).sub(L.benefit).sub(D.benefit);
 return{...d,L,D,p,A,C,value,benefit:letter.sub(value)};
}
function sfQSoftCommonInverseGuard(wt,x,mu){
 const m=wt.length,u=x.length,n=2*u-1,ny=m*(u+1),field=v=>[...v,sfdot(v,mu.slice(0,u-1)).neg().div(mu[u-1])],free=v=>{const l=field(v.slice(0,u-1)),d=field(v.slice(u-1,2*(u-1))),p=v.at(-1);return{A:l.map((z,i)=>z.add(d[i])),C:l.map((z,i)=>z.sub(d[i]).add(p))};},zero=x.map(()=>sfq(0)),ZY=Array.from({length:m},()=>zero),r0=wt.map(()=>sfq(0)),f0=v=>{const z=free(v);return sfQMeasureModelGeneralExpanded(wt,x,mu,z.A,z.C,ZY,r0).base;},M=sfQP6FormMatrix(n,f0),MI=sfQInv(M),basis=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>sfq(i===j?1:0))),fv=basis.map(f0),ys=Array.from({length:ny},(_,i)=>Array.from({length:ny},(_,j)=>sfq(i===j?1:0)));
 for(let i=0;i<ny-1;i++)ys.push(Array.from({length:ny},(_,j)=>sfq(j===i||j===i+1?1:0)));ys.push(Array.from({length:ny},(_,j)=>sfq(j%3-1)));
 let checked=0;for(const y of ys){const Z=Array.from({length:m},(_,j)=>y.slice(j*(u+1),j*(u+1)+u)),rho=Array.from({length:m},(_,j)=>y[j*(u+1)+u]),d=sfQSoftCommonInverse(wt,x,mu,Z,rho),letter=sfQMeasureModelGeneralExpanded(wt,x,mu,zero,zero,Z,rho).base,source=basis.map((v,i)=>{const z=free(v);return sfQMeasureModelGeneralExpanded(wt,x,mu,z.A,z.C,Z,rho).base.sub(fv[i]).sub(letter).div(2);}),explicit=letter.sub(sfdot(source,sfMatVec(MI,source))),atOpt=sfQMeasureModelGeneralExpanded(wt,x,mu,d.A,d.C,Z,rho).base;
 if(!d.value.eq(explicit)||!d.value.eq(atOpt))throw Error("soft common inverse formula against original matrix");checked++;}
 return{m,u,freeDimension:n,letterDimension:ny,exactInverseComparisons:checked,status:"PASS"};
}


function sfQCommonGaugeVerify(){
 const identities=[
 [[sfq(1),sfq(1,2)],[sfq(0),sfq(1,2),sfq(1)],[sfq(49,100),sfq(1,50),sfq(49,100)]],
 [[sfq(1,10),sfq(1,5),sfq(3,10)],[sfq(1,10),sfq(1,5),sfq(7,10)],[sfq(1,10),sfq(1,5),sfq(7,10)]],
 [[sfq(2),sfq(3,2),sfq(1,2)],[sfq(0),sfq(1,5),sfq(7,10),sfq(1)],[sfq(1,10),sfq(1,5),sfq(3,10),sfq(2,5)]]
 ].map(p=>sfQPCommonIdentityGuard(...p)),sources=[];
 for(const wU of [[sfq(3,5),sfq(4,5)],[sfq(1,3),sfq(2,3),sfq(2,3)],[sfq(1,2),sfq(1,2),sfq(1,2),sfq(1,2)]])for(const wt of [sfq(1,10),sfq(1),sfq(10)])sources.push(sfQPurePSourceGuard(wU,wt));
 const inverse=[
 [[sfq(1,100),sfq(1,100)],[sfq(1)],[sfq(1)]],
 [[sfq(1,2)],[sfq(0),sfq(1,2),sfq(1)],[sfq(49,100),sfq(1,50),sfq(49,100)]],
 [[sfq(1),sfq(1,2)],[sfq(1,10),sfq(1,5),sfq(7,10)],[sfq(1,10),sfq(1,5),sfq(7,10)]],
 [[sfq(1,10),sfq(1,5),sfq(3,10)],[sfq(0),sfq(1,5),sfq(7,10),sfq(1)],[sfq(1,10),sfq(1,5),sfq(3,10),sfq(2,5)]]
 ].map(p=>sfQSoftCommonInverseGuard(...p));
 const K=[sfq(1,2),sfq(1)],K1=spAdd(K,spC(1)),P=spSub(spScale(spProd(K,spSqCompat(K1),[sfq(1),sfq(3)]),sfq(4)),[sfq(0),sfq(0),sfq(2),sfq(2)]),expected=[sfq(9,2),sfq(57,2),sfq(57),sfq(44),sfq(12)];
 if(P.length!==expected.length||P.some((z,i)=>!z.eq(expected[i])||!z.gt(0)))throw Error("pure p half-operator polynomial");
 const wt=[sfq(1,100),sfq(1,100)],x=[sfq(1)],mu=x,r=[sfq(0),sfq(0)],both=sfQSoftCommonInverse(wt,x,mu,[[sfq(1)],[sfq(1)]],r),one=sfQSoftCommonInverse(wt,x,mu,[[sfq(1)],[sfq(0)]],r),two=sfQSoftCommonInverse(wt,x,mu,[[sfq(0)],[sfq(1)]],r),cross=both.benefit.sub(one.benefit).sub(two.benefit);
 if(!cross.gt(0))throw Error("independent-copy mean completion control did not fail");
 return{status:"PASS",identities,sources,inverse,operatorPolynomial:expected.map(z=>z.str()),independentMeanControlGap:cross.str()};
}
function spSqCompat(p){return spMul(p,p);}


function sfQSoftInverseSources(t,x,mu,JL,JD,b){
 const K=t.add(sfq(1,2)),f=x.map(s=>sfq(1).div(K.add(s))),F=sfdot(mu,f),alpha=sfq(1).div(sfq(2).add(F)),cp=t.sub(t.sq().mul(alpha).mul(F).mul(2)).sub(t.sq().div(K.sq().mul(2))),mL=x.map((s,i)=>s.add(sfq(1,2)).add(alpha.mul(s).mul(2)).mul(f[i]).mul(t).mul(2)),mD=x.map((s,i)=>s.add(sfq(1,2)).sub(alpha.mul(s).mul(2)).mul(f[i]).mul(t).mul(2)),cD=t.sq().mul(alpha).mul(2),cL=cD.sub(t.sq().sq().mul(alpha.sq()).mul(16).div(cp)),shift=t.sq().mul(alpha).mul(b).mul(4).div(cp),L=sfQSoftCommonCovariance(mL,f,JL.map((z,i)=>z.add(shift.mul(f[i]))),mu,cL),D=sfQSoftCommonCovariance(mD,f,JD,mu,cD),p=t.sq().mul(alpha).mul(sfdot(mu,L.optimizer.map((z,i)=>z.mul(f[i])))).mul(4).sub(b).div(cp),A=L.optimizer.map((z,i)=>z.add(D.optimizer[i])),C=L.optimizer.map((z,i)=>z.sub(D.optimizer[i]).add(p)),benefit=b.sq().div(cp).add(L.benefit).add(D.benefit);
 return{A,C,p,benefit,L,D,cp};
}


function sfQProjectedSource(wt,wU,A,C,Z,rho){
 const x=wU.map(w=>w.sq()),u=x.length,Pi=wU.map((w,i)=>wU.map((v,j)=>sfq(i===j?1:0).sub(w.mul(v)))),q=sfQMeasureModelGeneralExpanded(wt,x,x,A,C,Z,rho),D=x.map((s,i)=>x.map((_,j)=>i===j?s.mul(q.V[i]):sfq(0))),PD=sfQM(sfQM(Pi,D),Pi),R=PD.map((row,i)=>row.map((z,j)=>z.sub(wU[i].mul(wU[j]).mul(sfsum(q.P.map((p,l)=>p[i].mul(q.Q[l][j]))))))),norm=sfsum(R.flat().map(z=>z.sq()));
 if(!norm.eq(q.norm))throw Error("physical projected source norm");return{R,Pi,q};
}
function sfQPhysicalFreeSplitGuard(wt,wU){
 const t=sfdot(wt,wt),u=wU.length,x=wU.map(w=>w.sq()),A=x.map((_,i)=>sfq(i+1)),last=sfdot(x.slice(0,u-1),A.slice(0,u-1)).neg().div(x[u-1]);A[u-1]=last;const C=x.map((_,i)=>sfq(i-1,3)),Z=wt.map(()=>x.map(()=>sfq(0))),rho=wt.map(()=>sfq(0)),q=sfQProjectedSource(wt,wU,A,C,Z,rho),soft=wU.map((w,i)=>wU.map((v,j)=>w.mul(v).mul(A[i].add(C[j])))),base=q.q.base;
 if(!sfCUForm(wU,t,soft).eq(base))throw Error("physical soft map");
 let forward=0,transpose=0;for(let r=0;r<u;r++)for(let h=0;h<u;h++){const V=x.map((_,i)=>x.map((_,j)=>sfq(i===r&&j===h?1:0))),hard=sfHardProject(wU,V),actual=sfCUForm(wU,t,soft.map((row,i)=>row.map((v,j)=>v.add(hard[i][j])))).sub(base).sub(sfHardForm(wU,t,hard)),b=sfsum(hard.flatMap((row,i)=>row.map((v,j)=>v.mul(q.R[i][j])))).mul(2),bt=sfsum(hard.flatMap((row,i)=>row.map((v,j)=>v.mul(q.R[j][i])))).mul(2);if(actual.eq(b))forward++;if(actual.eq(bt))transpose++;}
 if(forward!==u*u&&transpose!==u*u)throw Error("physical source orientation");return{sigma:t.str(),u,sourceOrientation:forward===u*u?"as-written":"transpose",entryTests:u*u};
}


function sfQHardApply(wU,t,X){
 const x=wU.map(w=>w.sq()),K=t.add(sfq(1,2)),r=x.map(s=>s.div(K.add(s))),G=sfq(2).add(sfsum(r));
 return X.map((row,i)=>row.map((v,j)=>i===j?sfq(3,2).sub(r[i]).add(r[i].sq().div(G)).mul(v):v.div(2).add(r[i].mul(r[j]).div(G).mul(X[j][i]))));
}
function sfQHardTraceCalibration(wt,wU){
 const t=sfdot(wt,wt),K=t.add(sfq(1,2)),x=wU.map(w=>w.sq()),F=sfdot(x,x.map(s=>sfq(1).div(K.add(s)))),F2=sfdot(x,x.map(s=>sfq(1).div(K.add(s)).sq())),alpha=sfq(1).div(sfq(2).add(F)),JL=x.map(s=>t.mul(2).div(K.add(s)).mul(sfq(1).sub(alpha.mul(s).mul(2)).sub(alpha.mul(s).div(K.add(s))))),JD=x.map(()=>sfq(0)),b=t.sq().div(K).mul(K.add(1).mul(F).sub(1)).add(t.sq().mul(alpha).mul(F2.sub(F.sq()))),opt=sfQSoftInverseSources(t,x,x,JL,JD,b),zero=x.map(()=>sfq(0)),data=sfQProjectedSource(wt,wU,opt.A,opt.C,wt.map(()=>zero),wt.map(()=>sfq(0))),Pi=data.Pi,H=sfHardProject(wU,sfQHardApply(wU,t,Pi)),S=sfQAdd(H,data.R),pivot=sfsum(S.map((row,i)=>row[i])),expected=sfHardForm(wU,t,Pi).sub(opt.benefit);
 if(!pivot.eq(expected)||!pivot.gt(sfq(wU.length-1,72)))throw Error("hard trace calibration");
 return{t,x,Pi,S,pivot,opt,JL,b};
}
function sfQHardTraceLowerValue(wt,wU,Z,rho,constant=sfq(72)){
 const x=wU.map(w=>w.sq()),u=wU.length,soft=sfQSoftCommonInverse(wt,x,x,Z,rho);if(u===1)return{value:soft.value,soft:soft.value,residualNorm:sfq(0)};
 const c=sfQHardTraceCalibration(wt,wU),source=sfQProjectedSource(wt,wU,soft.A,soft.C,Z,rho),r0=sfsum(source.R.map((row,i)=>row[i])),res=source.R.map((row,i)=>row.map((v,j)=>v.sub(r0.div(c.pivot).mul(c.S[i][j])))),norm=sfsum(res.flat().map(v=>v.sq())),value=soft.value.sub(r0.sq().div(c.pivot)).sub(constant.mul(norm));
 if(!sfsum(res.map((row,i)=>row[i])).eq(0))throw Error("trace residual not trace-free");return{value,soft:soft.value,residualNorm:norm,pivot:c.pivot};
}
function sfQHardTraceLowerGuard(wt,wU){
 const u=wU.length,m=wt.length,t=sfdot(wt,wt),x=wU.map(w=>w.sq()),ns=2*u-1,nh=(u-1)**2,n=ns+nh,ny=m*(u+1),field=z=>[...z,sfdot(z,x.slice(0,u-1)).neg().div(x[u-1])],basis=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>sfq(i===j?1:0))),H=Array.from({length:u},(_,i)=>Array.from({length:u-1},(_,j)=>i===j?sfq(1):i===u-1?wU[j].neg().div(wU[u-1]):sfq(0)));
 function free(v){const l=field(v.slice(0,u-1)),d=field(v.slice(u-1,2*u-2)),p=v[2*u-2],A=l.map((z,i)=>z.add(d[i])),C=l.map((z,i)=>z.sub(d[i]).add(p)),vh=Array.from({length:u-1},(_,i)=>v.slice(ns+i*(u-1),ns+(i+1)*(u-1))),hard=sfQM(sfQM(H,vh),sfQT(H)),soft=wU.map((w,i)=>wU.map((z,j)=>w.mul(z).mul(A[i].add(C[j]))));return{A,C,hard,physical:sfQAdd(soft,hard)};}
 const Cfree=sfQP6FormMatrix(n,v=>sfCUForm(wU,t,free(v).physical)),CI=sfQInv(Cfree),zero=x.map(()=>sfq(0)),Z0=wt.map(()=>zero),rho0=wt.map(()=>sfq(0)),freevals=basis.map(v=>sfCUForm(wU,t,free(v).physical)),lower=[];
 for(let j=0;j<n;j++){const f=free(basis[j]),q=sfQProjectedSource(wt,wU,f.A,f.C,Z0,rho0),check=q.q.base.add(sfHardForm(wU,t,f.hard)).add(sfsum(f.hard.flatMap((row,i)=>row.map((v,k)=>v.mul(q.R[i][k])))).mul(2));if(!check.eq(freevals[j]))throw Error("free identity for complete trace guard");}
 function trueValue(y){const Z=Array.from({length:m},(_,j)=>y.slice(j*(u+1),j*(u+1)+u)),rho=Array.from({length:m},(_,j)=>y[j*(u+1)+u]),letter=sfQMeasureModelGeneralExpanded(wt,x,x,zero,zero,Z,rho).base,bb=basis.map((v,i)=>{const f=free(v),r=sfQProjectedSource(wt,wU,f.A,f.C,Z,rho),cross=sfsum(f.hard.flatMap((row,h)=>row.map((v,l)=>v.mul(r.R[h][l])))).mul(2),all=r.q.base.add(sfHardForm(wU,t,f.hard)).add(cross);return all.sub(freevals[i]).sub(letter).div(2);});return letter.sub(sfdot(bb,sfMatVec(CI,bb)));}
 const G=sfQP6FormMatrix(ny,trueValue),L=sfQP6FormMatrix(ny,y=>{const Z=Array.from({length:m},(_,j)=>y.slice(j*(u+1),j*(u+1)+u)),rho=Array.from({length:m},(_,j)=>y[j*(u+1)+u]);return sfQHardTraceLowerValue(wt,wU,Z,rho).value;}),I=sfQInertia(sfQSub(G,L));
 if(I.negative||(u===2&&(I.positive||I.zero!==ny)))throw Error("hard trace lower form against exact free Schur");
 return{sigma:t.str(),u,m,freeDimension:n,letterDimension:ny,freeIdentityChecks:n,gap:I,lower:sfQInertia(L)};
}


function sfQSoftHardShellGuard(wt,wU){
 const u=wU.length,t=sfdot(wt,wt),x=wU.map(w=>w.sq()),ns=2*u-1,nh=(u-1)**2,field=v=>[...v,sfdot(v,x.slice(0,u-1)).neg().div(x[u-1])],softBasis=Array.from({length:ns},(_,j)=>{const v=Array.from({length:ns},(_,i)=>sfq(i===j?1:0)),L=field(v.slice(0,u-1)),D=field(v.slice(u-1,2*u-2)),p=v.at(-1),A=L.map((z,i)=>z.add(D[i])),C=L.map((z,i)=>z.sub(D[i]).add(p)),physical=wU.map((w,i)=>wU.map((v,j)=>w.mul(v).mul(A[i].add(C[j]))));return{A,C,physical};}),H=Array.from({length:u},(_,i)=>Array.from({length:u-1},(_,j)=>i===j?sfq(1):i===u-1?wU[j].neg().div(wU[u-1]):sfq(0))),hardBasis=Array.from({length:nh},(_,j)=>{const E=Array.from({length:u-1},(_,r)=>Array.from({length:u-1},(_,c)=>sfq(r*(u-1)+c===j?1:0)));return sfQM(sfQM(H,E),sfQT(H));}),zero=x.map(()=>sfq(0)),R=softBasis.map(v=>sfQProjectedSource(wt,wU,v.A,v.C,wt.map(()=>zero),wt.map(()=>sfq(0))).R),inner=(A,B)=>sfsum(A.flatMap((row,i)=>row.map((v,j)=>v.mul(B[i][j])))),sumBasis=(basis,v)=>basis[0].map((row,i)=>row.map((_,j)=>sfsum(basis.map((B,k)=>B[i][j].mul(v[k]))))),A=sfQP6FormMatrix(ns,v=>sfCUForm(wU,t,sumBasis(softBasis.map(v=>v.physical),v))),Ainv=sfQInv(A),T=hardBasis.map(X=>R.map(Y=>inner(X,Y))),gram=hardBasis.map(X=>hardBasis.map(Y=>inner(X,Y))),hard=sfQP6FormMatrix(nh,v=>sfHardForm(wU,t,sumBasis(hardBasis,v))),loss=sfQM(T,sfQM(Ainv,sfQT(T))),relative=sfQInertia(sfQSub(sfQScale(gram,sfq(29,72)),loss)),completed=sfQInertia(sfQSub(sfQSub(hard,loss),sfQScale(gram,sfq(1,72))));
 if(relative.negative||completed.negative||relative.zero||completed.zero)throw Error("soft hard static shell");
 return{sigma:t.str(),u,softDimension:ns,hardDimension:nh,relative,completed};
}


function sfQKernelNorm(K,x){return sfsum(K.flatMap((row,i)=>row.map((v,j)=>x[i].mul(x[j]).mul(v.sq()))));}
function sfQKernelCenter(K,x){const rows=K.map(row=>sfdot(x,row)),cols=K[0].map((_,j)=>sfdot(x,K.map(row=>row[j]))),mean=sfdot(x,rows);return K.map((row,i)=>row.map((v,j)=>v.sub(rows[i]).sub(cols[j]).add(mean)));}
function sfQSourceKernel(wt,x,A,C,Z,rho){
 const d=sfQMeasureModelGeneralExpanded(wt,x,x,A,C,Z,rho),E=sfdot(x,x.map((s,i)=>s.mul(d.V[i]))),K=x.map((s,i)=>x.map((ss,j)=>sfq(i===j?1:0).mul(d.V[i]).sub(s.mul(d.V[i])).sub(ss.mul(d.V[j])).add(E).sub(sfsum(d.P.map((p,l)=>p[i].mul(d.Q[l][j]))))));
 if(!sfQKernelNorm(K,x).eq(d.norm))throw Error("weighted physical kernel norm");return{K,d};
}
function sfQWeightedTraceCalibration(wt,x){
 const t=sfdot(wt,wt),K=t.add(sfq(1,2)),f=x.map(s=>sfq(1).div(K.add(s))),F=sfdot(x,f),F2=sfdot(x,f.map(s=>s.sq())),alpha=sfq(1).div(sfq(2).add(F)),JL=x.map((s,i)=>t.mul(2).mul(f[i]).mul(sfq(1).sub(alpha.mul(s).mul(2)).sub(alpha.mul(s).mul(f[i])))),JD=x.map(()=>sfq(0)),b=t.sq().div(K).mul(K.add(1).mul(F).sub(1)).add(t.sq().mul(alpha).mul(F2.sub(F.sq()))),opt=sfQSoftInverseSources(t,x,x,JL,JD,b),zero=x.map(()=>sfq(0)),src=sfQSourceKernel(wt,x,opt.A,opt.C,wt.map(()=>zero),wt.map(()=>sfq(0))),r=x.map((s,i)=>s.mul(f[i])),G=sfq(2).add(F),Pi=x.map((s,i)=>x.map((ss,j)=>sfq(i===j?1:0).div(s).sub(1))),HPi=Pi.map((row,i)=>row.map((v,j)=>i===j?sfq(3,2).sub(r[i]).add(r[i].sq().div(G)).mul(v):v.div(2).add(r[i].mul(r[j]).div(G).mul(Pi[j][i])))),Hproj=sfQKernelCenter(HPi,x),S=sfQAdd(Hproj,src.K),pivot=sfdot(x,S.map((row,i)=>row[i])),hard=sfsum(Pi.flatMap((row,i)=>row.map((v,j)=>v.mul(HPi[i][j]).mul(x[i]).mul(x[j]))));
 if(!pivot.eq(hard.sub(opt.benefit))||!pivot.gt(sfq(x.length-1,72)))throw Error("weighted trace pivot");return{S,pivot};
}
function sfQWeightedTraceLower(wt,x,Z,rho){
 const soft=sfQSoftCommonInverse(wt,x,x,Z,rho);if(x.length===1)return{value:soft.value};
 const trace=sfQWeightedTraceCalibration(wt,x),src=sfQSourceKernel(wt,x,soft.A,soft.C,Z,rho),r0=sfdot(x,src.K.map((row,i)=>row[i])),res=src.K.map((row,i)=>row.map((v,j)=>v.sub(r0.div(trace.pivot).mul(trace.S[i][j]))));
 if(!sfdot(x,res.map((row,i)=>row[i])).eq(0))throw Error("weighted trace residual");const norm=sfQKernelNorm(res,x),value=soft.value.sub(r0.sq().div(trace.pivot)).sub(norm.mul(72));return{value,soft:soft.value,traceGain:r0.sq().div(trace.pivot),traceFreeSourceNorm:norm};
}
function sfQTraceCoarseBoundControl(){
 const wt=[sfq(1,100),sfq(1,10),sfq(1)],x=[sfq(1,10000),sfq(9,10000),sfq(9,1000),sfq(9,100),sfq(3,10),sfq(3,5)],zero=x.map(()=>sfq(0)),Z=[zero,zero,x.map(()=>sfq(1))],rho=[sfq(0),sfq(0),sfq(1)],r=sfQWeightedTraceLower(wt,x,Z,rho);
 if(!sfq(-3,8).gt(r.value))throw Error("coarse trace relaxation did not fail");
 return{status:"EXPECTED_FAILURE_CERTIFIED",sigma:sfdot(wt,wt).str(),squaredFreeWeights:x.map(z=>z.str()),lowerValue:r.value.str(),approx:r.value.num(),strictBound:"< -3/8",scope:"negative sufficient lower bound only; no NPT witness"};
}


function sfQCommonGaugeAllVerify(){
 const basic=sfQCommonGaugeVerify(),profiles=[
 [[sfq(1,10),sfq(1)],[sfq(3,5),sfq(4,5)]],
 [[sfq(1),sfq(1,2)],[sfq(1,3),sfq(2,3),sfq(2,3)]],
 [[sfq(2),sfq(3)],[sfq(1,2),sfq(1,2),sfq(1,2),sfq(1,2)]]
 ],physical=profiles.map(p=>sfQPhysicalFreeSplitGuard(...p)),shell=profiles.map(p=>sfQSoftHardShellGuard(...p)),trace=[
 [[sfq(1,10),sfq(1)],[sfq(3,5),sfq(4,5)]],
 [[sfq(1),sfq(1,2)],[sfq(1,3),sfq(2,3),sfq(2,3)]],
 [[sfq(1,10),sfq(1,5),sfq(3,10)],[sfq(1,3),sfq(2,3),sfq(2,3)]]
 ].map(p=>sfQHardTraceLowerGuard(...p));
 let kernelIdentities=0;for(const [wt,wU]of profiles){const r=sfQHardTraceCalibration(wt,wU),k=sfQWeightedTraceCalibration(wt,wU.map(w=>w.sq()));if(!r.pivot.eq(k.pivot))throw Error("weighted trace pivot");for(let i=0;i<wU.length;i++)for(let j=0;j<wU.length;j++){if(!r.S[i][j].eq(wU[i].mul(wU[j]).mul(k.S[i][j])))throw Error("weighted trace action");kernelIdentities++;}}
 const coarse=sfQTraceCoarseBoundControl();return{status:"PASS",basic,physical,shell,trace,kernelIdentities,coarseControl:{...coarse,lowerValue:undefined}};
}


function sfQReviewA510RelativeDiagonalGuard(x,V){
 const E=y=>sfdot(x,y),c=sfq(5,18),b=sfq(13,72),d=sfq(11,2592),B=E(x.map((s,i)=>s.add(sfq(1,4)).mul(V[i].sq()))).add(E(V).sq().div(4)),N=E(x.map((s,i)=>s.mul(sfq(1).sub(s.mul(2))).mul(V[i].sq()))).add(E(x.map((s,i)=>s.mul(V[i]))).sq()),Z=x.map((s,i)=>s.sub(b).mul(V[i])),variance=E(Z.map(v=>v.sq())).sub(E(Z).sq()),rhs=variance.mul(2).add(E(x.map((s,i)=>s.sub(b.mul(2)).mul(V[i]))).sq()).add(d.mul(E(V.map(v=>v.sq())).add(E(V).sq())));
 if(!c.mul(B).sub(N).eq(rhs)||!rhs.ge(0))throw Error("relative diagonal rational identity");return{B,N,rhs};
}
function sfQReviewA510RelaxGuard(w,wt,seed=1){
 const u=w.length,t=wt.sq(),x=w.map(z=>z.sq()),K=t.add(sfq(1,2)),k=t.add(1).div(2),ww=[wt,...w],W=sfInvE(ww,k),Z=x.map((_,i)=>sfq((i+1)*seed-2,i+2)),rho=sfq(seed-2,3),zeta=w.map((z,i)=>z.mul(Z[i])),Y=[rho,...zeta],WY=sfMatVec(W,Y),Y0=wt.mul(rho).add(sfdot(w,zeta)),m=t.mul(2).add(1).mul(rho).add(wt.mul(sfdot(w,zeta))).div(t.add(1)),Aletter=rho.sq().mul(2).add(sfdot(zeta,zeta)).sub(m.sq().div(2)),cr=w.map((z,i)=>Y0.div(2).add(z.mul(zeta[i]))),H=w.map((wr,r)=>w.map((wh,h)=>sfq(r===h?1:0).mul(rho.sub(m)).add(m.mul(wr).mul(wh).div(K.mul(2))).sub(cr[r].mul(W[0][r+1]).mul(wh)).add(wt.mul(cr[r]).mul(W[h+1][r+1])).sub(WY[0].mul(wr).mul(wh).div(2)).add(wt.mul(wr).mul(WY[h+1]).div(2)))),A=x.map((_,i)=>sfq(i*seed+1,i+3));A[u-1]=sfdot(x.slice(0,-1),A.slice(0,-1)).neg().div(x[u-1]);const C=x.map((_,i)=>sfq(seed-i,5)),Omega=w.map((wr,r)=>w.map((wh,h)=>wr.mul(wh).mul(A[r].add(C[h])))),data=sfQProjectedSource([wt],w,A,C,[Z],[rho]),base=data.q.base,noteBase=X=>sfCUForm(w,t,X).add(sfsum(H.flatMap((row,i)=>row.map((v,j)=>v.mul(X[i][j]))))).add(Aletter),forward=base.eq(noteBase(Omega)),backward=base.eq(noteBase(sfQT(Omega)));
 if(!forward&&!backward)throw Error("a510 q-eliminated normalization");
 const X=forward?Omega:sfQT(Omega),av=sfMatVec(X,w),cv=sfMatVec(sfQT(X),w),p=sfdot(w,cv),aa=av.map((z,i)=>z.sub(w[i].mul(sfdot(w,av)))),Ei=sfInvE(w,K),G=w.map((_,i)=>w.map((_,j)=>Ei[i][j].mul(aa[j]).mul(w[j]).mul(2).add(cv[i].mul(w[i]).mul(Ei[i][j]).mul(2)).add(p.mul(Ei[i][j].sub(sfq(i===j?1:0).div(K)))))),total=sfHardProject(w,G.map((row,i)=>row.map((v,j)=>v.mul(t).add(H[i][j])))),R=data.R,Rt=sfQT(R),matches=total.every((row,i)=>row.every((v,j)=>v.eq(R[i][j].mul(2)))),matchesTranspose=total.every((row,i)=>row.every((v,j)=>v.eq(Rt[i][j].mul(2))));
 if(!matches&&!matchesTranspose)throw Error("a510 full hard source factor two");
 const aa0=sfQSub(sfQDiag(w.map(()=>sfq(1))),sfQScale(Ei,t)),zv=sfMatVec(Ei,w),Mac=w.map((_,i)=>w.map((_,j)=>zv[i].mul(w[j]).div(2).add(i===j?w[i].mul(zv[i]):sfq(0)))),Mcc=aa0.map((row,i)=>row.map((v,j)=>v.sub(t.mul(w[i]).mul(w[j]).div(K.sq().mul(2))))),soft=t.mul(sfdot(aa,sfMatVec(aa0,aa)).add(sfdot(aa,sfMatVec(Mac,cv)).mul(2)).add(sfdot(cv,sfMatVec(Mcc,cv))));
 if(!soft.eq(sfCUForm(w,t,X)))throw Error("a510 soft F2 identity");
 const lambda=w.map((_,i)=>w.map((_,j)=>i===j?sfq(3,2).sub(x[i].div(K.add(x[i]))):sfq(5,12))),loss=sfsum(total.flatMap((row,i)=>row.map((v,j)=>v.sq().div(lambda[i][j].mul(4))))),ours=sfsum(R.flatMap((row,i)=>row.map((v,j)=>v.sq().div(lambda[i][j]))));
 if(!loss.eq(ours))throw Error("a510 entrywise relaxation identity");
 const triangle=sfQTriangleRemainder([wt],x,x,A,C,[Z],[rho]),diag=sfQReviewA510RelativeDiagonalGuard(x,data.q.V),norm=sfsum(R.flat().map(v=>v.sq())),T=triangle.B.add(triangle.R.mul(8)),Psi=base.sub(loss),comparison=base.sub(T.mul(sfq(29,30)));
 if(!triangle.base.eq(base)||!diag.B.eq(triangle.B)||!T.mul(sfq(29,72)).ge(norm)||!norm.mul(sfq(12,5)).ge(loss)||!triangle.remainder.gt(0)||!Psi.ge(comparison)||!comparison.gt(base.div(30)))throw Error("a510 Conjecture R comparison guard");
 return{u,sigma:t.str(),seed,baseConvention:forward?"as-written":"transpose",sourceConvention:matches?"as-written":"transpose",entrywiseIdentity:true,positiveEnergy:true,PsiPositive:true,relativeMarginAboveBaseOverThirty:comparison.sub(base.div(30)).gt(0)};
}
function sfQReviewA510RelaxAllVerify(){
 const profiles=[[sfq(3,5),sfq(4,5)],[sfq(2,7),sfq(3,7),sfq(6,7)],[sfq(1,5),sfq(2,5),sfq(2,5),sfq(4,5)],[sfq(2,5),sfq(2,5),sfq(2,5),sfq(2,5),sfq(3,5)]],cases=[];
 for(const w of profiles)for(const wt of [sfq(1,10),sfq(1),sfq(3,2)])cases.push(sfQReviewA510RelaxGuard(w,wt,cases.length%3+1));
 return{status:"PASS",cases,scope:"Independent rational checks of the new Lemma R1 normalization, full projected source, entrywise relaxation, and its strict lower bound from the already proved one-letter energy theorem."};
}

