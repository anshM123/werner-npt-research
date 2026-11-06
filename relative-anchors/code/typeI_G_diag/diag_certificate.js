function sfP6BuildCertificate(ax,gi){

const MPBASE=256,MPSQ=65536;let mpOps=0,mpMaxTerms=0,mpStart=Date.now(),mpDeadline=mpStart+45000;
const mpGcd=(a,b)=>{a=a<0n?-a:a;b=b<0n?-b:b;while(b){const c=a%b;a=b;b=c;}return a;};
const mpConst=n=>BigInt(n)===0n?new Map():new Map([[0,BigInt(n)]]),mpVar=i=>new Map([[i===0?MPSQ:i===1?MPBASE:1,1n]]);
const mpTerms=p=>{mpMaxTerms=Math.max(mpMaxTerms,p.size);if((++mpOps&255)===0&&Date.now()>mpDeadline)throw Error("MP time guard");return p;};
const mpAdd=(p,q)=>{const r=new Map(p);for(const[k,v]of q){const n=(r.get(k)||0n)+v;if(n)r.set(k,n);else r.delete(k);}return mpTerms(r);};
const mpScale=(p,c)=>c===0n?new Map():new Map([...p].map(([k,v])=>[k,v*c])),mpNeg=p=>mpScale(p,-1n),mpSub=(p,q)=>mpAdd(p,mpNeg(q));
function mpMul(p,q){if(!p.size||!q.size)return new Map();const dp=mpDegree(p),dq=mpDegree(q);if(dp.some((v,i)=>v+dq[i]>=MPBASE))throw Error("MP exponent encoding guard");const r=new Map();for(const[k,v]of p)for(const[l,w]of q){const z=k+l,c=(r.get(z)||0n)+v*w;if(c)r.set(z,c);else r.delete(z);}return mpTerms(r);}
const mpContent=p=>{let c=0n;for(const v of p.values())c=mpGcd(c,v);return c;},mpDivScalar=(p,d)=>new Map([...p].map(([k,v])=>{if(v%d)throw Error("MP scalar nonexact");return[k,v/d];}));
const mpLead=p=>{let k=-1;for(const l of p.keys())if(l>k)k=l;return[k,p.get(k)];};
const mpExp=k=>[Math.floor(k/MPSQ),Math.floor(k/MPBASE)%MPBASE,k%MPBASE];
const mpDegree=p=>{const d=[0,0,0];for(const k of p.keys()){const e=mpExp(k);for(let i=0;i<3;i++)d[i]=Math.max(d[i],e[i]);}return d;};

function mpTryDiv(p,q){if(!q.size)throw Error("MP zero divide");if(!p.size)return new Map();const[d,lc]=mpLead(q),de=mpExp(d),sumq=[...q.values()].reduce((a,b)=>a+b,0n),sump=[...p.values()].reduce((a,b)=>a+b,0n);if(sumq&&sump%sumq)return null;
 const r=new Map(p),out=new Map(),heap=[...p.keys()];heap.sort((a,b)=>b-a);
 const push=v=>{let i=heap.length;heap.push(v);while(i){const j=(i-1)>>1;if(heap[j]>=v)break;heap[i]=heap[j];i=j;}heap[i]=v;};
 const pop=()=>{const ret=heap[0],last=heap.pop();if(heap.length){let i=0;while(i*2+1<heap.length){let j=i*2+1;if(j+1<heap.length&&heap[j+1]>heap[j])j++;if(heap[j]<=last)break;heap[i]=heap[j];i=j;}heap[i]=last;}return ret;};
 let steps=0;while(heap.length){const key=pop(),val=r.get(key);if(!val)continue;const ke=mpExp(key);if(ke.some((x,i)=>x<de[i])||val%lc)return null;const h=key-d,c=val/lc;out.set(h,c);
 for(const[l,w]of q){const z=l+h,old=r.get(z)||0n,n=old-c*w;if(n){r.set(z,n);if(!old)push(z);}else r.delete(z);}
 if((++steps&4095)===0&&Date.now()>mpDeadline)throw Error("MP division time guard");}
 return mpTerms(out);
}

let mpAtoms=[],mpAtomKeys=new Map();
function mpKey(p){return [...p].sort((a,b)=>a[0]-b[0]).map(([k,v])=>k+":"+v).join(",");}
function mpRegister(p,label=""){const key=mpKey(p);if(mpAtomKeys.has(key))return mpAtomKeys.get(key);const id=mpAtoms.length;mpAtoms.push({p:new Map(p),label:label||"factor"+id,positive:[...p.values()].every(v=>v>=0n),degree:mpDegree(p)});mpAtomKeys.set(key,id);return id;}
const mpDenUnion=(x,y)=>{const d=new Map(x);for(const[k,v]of y)d.set(k,Math.max(v,d.get(k)||0));return d;},mpDenSum=(x,y)=>{const d=new Map(x);for(const[k,v]of y)d.set(k,v+(d.get(k)||0));return d;};
function mpApplyFactors(p,den,minus=new Map()){for(const[id,n]of den){const k=n-(minus.get(id)||0);for(let i=0;i<k;i++)p=mpMul(p,mpAtoms[id].p);}return p;}
class SFMPR{
 constructor(n,d=1n,atoms=new Map(),normalize=true){if(d===0n)throw Error("MP scalar zero denom");if(!n.size){this.n=new Map();this.d=1n;this.f=new Map();return;}n=new Map(n);atoms=new Map(atoms);if(d<0n){d=-d;n=mpNeg(n);}const c=mpGcd(mpContent(n),d);if(c>1n){n=mpDivScalar(n,c);d/=c;}if(normalize)for(const[id,num]of atoms){let left=num;while(left){const q=mpTryDiv(n,mpAtoms[id].p);if(!q)break;n=q;left--;}if(left)atoms.set(id,left);else atoms.delete(id);}this.n=n;this.d=d;this.f=atoms;}
 add(y){y=sfq(y);const den=mpDenUnion(this.f,y.f),g=mpGcd(this.d,y.d),d=this.d/g*y.d,n=mpAdd(mpScale(mpApplyFactors(this.n,den,this.f),d/this.d),mpScale(mpApplyFactors(y.n,den,y.f),d/y.d));return new SFMPR(n,d,den);}
 sub(y){return this.add(sfq(y).neg());}neg(){return new SFMPR(mpNeg(this.n),this.d,this.f,false);}
 mul(y){y=sfq(y);return new SFMPR(mpMul(this.n,y.n),this.d*y.d,mpDenSum(this.f,y.f));}
 div(y){y=sfq(y);if(!y.n.size)throw Error("MP divide zero");let p=new Map(y.n),c=mpContent(p);if(mpLead(p)[1]<0n)c=-c;p=mpDivScalar(p,c);const f=new Map();for(let id=0;id<mpAtoms.length&&p.size&&!(p.size===1&&p.has(0));id++){let num=0;while(true){const q=mpTryDiv(p,mpAtoms[id].p);if(!q)break;p=q;num++;}if(num)f.set(id,num);}if(!(p.size===1&&p.get(0)===1n)){const id=mpRegister(p);f.set(id,(f.get(id)||0)+1);}
 return this.mul(new SFMPR(mpApplyFactors(mpConst(y.d),y.f),c,f));}
 sq(){return this.mul(this);}
 eq(y){y=this.sub(y);return !y.n.size;}
}
const sfq=(a,b=1)=>a instanceof SFMPR?a:new SFMPR(mpConst(a),BigInt(b)),sfsum=a=>a.reduce((t,x)=>t.add(x),sfq(0)),sfdot=(a,b)=>sfsum(a.map((x,i)=>x.mul(b[i])));
function mpInitWeights(){const a=new SFMPR(mpVar(0)),env=new SFMPR(mpVar(1)),s=new SFMPR(mpVar(2)),t=a.add(env),K=t.add(sfq(1,2)),k=K.sub(a.div(2)),delta=k.add(a);for(const[name,x]of Object.entries({a,t,K,K1:K.add(1),k,k1:k.add(1),delta,t1:t.add(1),d:K.add(s),e:k.add(s)})){const c=mpContent(x.n);mpRegister(mpDivScalar(x.n,c),name);}return{a,env,s,t,K,k,delta};}
function mpTriangle(k,c){const lo=sfq(1).div(k.add(1)),hi=sfq(1).div(k),D=c.mul(k).mul(k.add(1)).mul(2).add(k.mul(2)).add(1);return[[lo,sfq(1).div(c.add(lo))],[hi,sfq(1).div(c.add(hi))],[c.mul(k.mul(2).add(1)).add(2).div(D),k.mul(k.add(1)).mul(2).div(D)]];}
function mpSerialize(p){return[...p].map(([k,v])=>[k,v.toString()]);}
function mpDeserialize(p){return new Map(p.map(([k,v])=>[k,BigInt(v)]));}

function sfQP6(t,a,s,F,alpha,B,gamma,H=sfq(3,4),centerMultiplier=sfq(1)){
 const K=t.add(sfq(1,2)),k=K.sub(a.div(2)),delta=k.add(a),theta=a.div(t),env=t.sub(a),u=t.mul(2).sub(a),v=k.sub(sfq(1,2)),D0=sfq(2).add(a.div(delta)),C0=sfq(2).add(a.add(1).div(delta)),d=K.add(s),e=k.add(s),f=sfq(1).div(d),g=sfq(1).div(e),om=s.add(sfq(1,4)),VF=sfq(1).div(K.sq().mul(K.add(1).sq()).mul(4)),VG=sfq(1).div(k.sq().mul(k.add(1).sq()).mul(4)),sqrtVF=sfq(1).div(K.mul(K.add(1)).mul(2)),sqrtVG=sfq(1).div(k.mul(k.add(1)).mul(2));
 const P=Array.from({length:6},()=>Array.from({length:6},()=>sfq(0))),add=(i,j,z)=>{P[i][j]=P[i][j].add(z);if(i!==j)P[j][i]=P[j][i].add(z);},outer=(co,x)=>{for(let i=0;i<6;i++)for(let j=0;j<6;j++)P[i][j]=P[i][j].add(co.mul(x[i]).mul(x[j]));};
 P[0][0]=s.add(sfq(1,2)).add(alpha.mul(s).mul(2)).mul(2).div(d).sub(t.mul(om).mul(4).div(d.sq()));
 P[1][1]=s.add(sfq(1,2)).sub(alpha.mul(s).mul(2)).mul(2).div(d);
 P[2][2]=sfq(1).sub(a.mul(om).div(e.sq().mul(4))).mul(t).div(a);
 add(0,2,gamma.mul(C0).mul(s).add(s.add(sfq(1,2))).sub(t.mul(om).mul(2).div(d)).div(e.mul(2)));
 add(1,2,gamma.mul(C0).mul(s).sub(s.add(sfq(1,2))).div(e.mul(2)));
 P[3][3]=sfq(1).sub(t.mul(alpha).mul(f).mul(2)).sub(t.div(K.sq().mul(2)));
 add(3,4,sfq(1).sub(sfq(1).div(delta.mul(2))).add(gamma.mul(a).div(2).mul(D0.mul(g).add(sfq(1).div(delta.sq())))).sub(t.div(k)).div(a.mul(2)));
 P[4][4]=sfq(2).sub(K.sq().div(k.sq().mul(2))).mul(t).div(a.sq());
 add(4,5,u.div(a).mul(sfq(2).sub(K.sq().div(k.sq().mul(2)))).sub(K.div(k.sq().mul(4))).mul(t).div(a));
 P[5][5]=sfq(1).add(env.mul(t.mul(3).sub(a)).mul(2).div(a)).mul(t).div(a);
 add(3,5,t.div(K.mul(a).mul(2)).sub(env.div(a.mul(2)).mul(a.add(1).div(delta).sub(gamma.mul(a).mul(sfq(1).div(delta.sq()).add(D0.div(delta)))))));
 const gammaMax=sfq(1).div(D0.add(sfq(1).div(k.add(1)))),cmax=a.mul(gammaMax).div(2).mul(a.div(k.mul(delta)).add(sfq(1).div(k.mul(k.add(1))))),epsEll=cmax.mul(sqrtVG),cEll=sfq(1).sub(theta).div(4),cE=cEll.mul(sfq(1).add(epsEll)).mul(t).div(a);
 function source(z){
   const ff=sfq(1).div(K.add(z)),gg=sfq(1).div(k.add(z)),oo=z.add(sfq(1,4));
   const WA=[[t.mul(alpha).mul(ff).mul(-4),gamma.mul(t.mul(2).add(a.mul(2)).add(1)).mul(gg).div(delta.mul(2)),gamma.mul(env).mul(gg).div(delta.mul(-2))],[sfq(0),gamma.mul(gg).div(delta.mul(2)),gamma.mul(env).mul(gg).div(delta.mul(2))],[gamma.mul(C0).mul(v).mul(gg).div(-2),sfq(0),sfq(0)]];
   const WB=[[t.sq().mul(2).div(K).mul(ff).mul(oo.mul(ff).add(F.div(4))),t.div(k.mul(-2)).mul(ff).mul(oo.mul(gg).add(B.div(4))),sfq(0)],[sfq(0),sfq(0),sfq(0)],[t.sq().div(K.mul(2)).mul(gg).mul(oo.mul(ff).add(F.div(4))),t.div(k.mul(-8)).mul(gg).mul(oo.mul(gg).add(B.div(4))),sfq(0)]];
   const ae=sfq(1).add(a.mul(gg).div(2)),be=sfq(1).sub(a.mul(gamma).div(delta)).add(a.mul(gamma).mul(gg)),mc=[sfq(0),sfq(1).div(delta),t.mul(3).add(1).div(delta)];
   return WA.map((row,i)=>row.map((x,j)=>x.sub(WB[i][j]).sub(i===2?cE.mul(ae).mul(be).mul(mc[j]):sfq(0))));
 }
 const W=source(s),W0=source(sfq(0));for(let i=0;i<3;i++)for(let j=0;j<3;j++)add(i,3+j,W[i][j].sub(W0[i][j].mul(centerMultiplier)).add(i===0&&j===0?t.div(t.add(1).sq().mul(t.add(1)).mul(4)):sfq(0)));
 const bs=[sfq(0),sfq(0),sfq(0),t.div(K).mul(f),g.div(k.mul(-4)),sfq(0)],bbs=[sfq(0),sfq(0),sfq(0),t.div(K).mul(F),B.div(k.mul(-4)),sfq(0)];outer(t.mul(om).neg(),bs);outer(t.div(-4),bbs);
 const zLoss=a.mul(VG).div(16).add(a.div(t).mul(K.add(1).sq()).mul(VG).mul(gamma.mul(v).mul(sqrtVG).mul(2).add(t.mul(sqrtVF)).sq())).add(gamma.sq().mul(a).div(t).mul(v.sq()).mul(K.add(1)).mul(VG.sq()).mul(4));
 P[0][0]=P[0][0].sub(t.mul(sfq(1).sub(alpha.mul(2))).div(K.sq().mul(4)).add(sfq(1,16)).mul(f.sq()));P[1][1]=P[1][1].sub(f.div(16));P[2][2]=P[2][2].sub(zLoss.mul(t).div(a));
 outer(t.mul(alpha.sq()).mul(VF).mul(-64),[s.mul(f),sfq(0),sfq(0),t.mul(f).div(-2),sfq(0),sfq(0)]);
 outer(t.mul(alpha.sq()).mul(VF).mul(-32),[sfq(0),s.mul(f),sfq(0),sfq(0),sfq(0),sfq(0)]);
 outer(t.mul(gamma.sq()).mul(VG).mul(-4),[sfq(0),sfq(0),s.mul(g),sfq(0),g.div(2),sfq(0)]);
 const ae=sfq(1).add(a.mul(g).div(2)),be=sfq(1).sub(a.mul(gamma).div(delta)).add(a.mul(gamma).mul(g)),mc=[sfq(0),sfq(1).div(delta),t.mul(3).add(1).div(delta)];
 outer(cE.neg(),[sfq(0),sfq(0),ae,sfq(0),be.div(delta),be.mul(t.mul(3).add(1)).div(delta)]);
 for(let j=0;j<3;j++)add(2,3+j,cE.mul(ae).mul(be).mul(mc[j]));
 P[2][2]=P[2][2].sub(cEll.mul(epsEll.add(epsEll.sq())).mul(t).div(a));
 const co=H.mul(sfq(1).sub(theta));P[2][2]=P[2][2].sub(co.mul(t).div(a));
 outer(co.mul(t).div(a.sq()).neg(),[sfq(0),sfq(0),sfq(0),sfq(0),sfq(1),u]);P[5][5]=P[5][5].sub(co.mul(t.sq()).div(a.mul(t.add(1))));
 return{P,K,k,delta,epsEll};
}

function mpLoadPoint(data){mpAtoms=[];mpAtomKeys=new Map();for(const x of data.atomData)mpRegister(mpDeserialize(x.p),x.label);return data.data.map(row=>row.map(x=>new SFMPR(mpDeserialize(x.n),BigInt(x.d),new Map(x.f))));}
function mpDet(M){const n=M.length,memo=new Map();function rec(mask){if(mask===(1<<n)-1)return sfq(1);if(memo.has(mask))return memo.get(mask);let row=0;for(let j=0;j<n;j++)if(mask>>j&1)row++;let z=sfq(0),sgn=1;for(let j=0;j<n;j++)if(!(mask>>j&1)){if(M[row][j].n.size){let v=M[row][j].mul(rec(mask|1<<j));if(sgn<0)v=v.neg();z=z.add(v);}sgn=-sgn;}memo.set(mask,z);return z;}return rec(0);}
function mpStats(p){let positive=0,negative=0,strictEdge=false;for(const[k,v]of p){if(v>0n)positive++;else if(v<0n)negative++;const e=mpExp(k);if(!e[1]&&!e[2]&&v>0n)strictEdge=true;}return{terms:p.size,positive,negative,degree:mpDegree(p),strictEdge};}
function mpRFSerialize(x){return{n:mpSerialize(x.n),d:x.d.toString(),f:[...x.f]};}

function mpP6Build(ax,gi,wrongCenter=false){
 const {a,env,s,t,K,k,delta}=mpInitWeights(),A=mpTriangle(K,sfq(2)),G=mpTriangle(k,sfq(2).add(a.div(delta))),P=sfQP6(t,a,s,A[ax][0],A[ax][1],G[gi][0],G[gi][1],sfq(3,4),sfq(wrongCenter?2:1)).P,S=[K.add(s),K.add(s),k.add(s),K.add(s),k.add(s),k.add(s)],Q=P.map((row,i)=>row.map((v,j)=>v.mul(S[i]).mul(S[j]))),ord=[5,1,2,4,0,3],minors=[];
 if(Q.some(row=>row.some(x=>[...x.f].some(([id])=>!mpAtoms[id].positive||mpAtoms[id].degree[2]))))throw Error("nonpositive or point-dependent denominator atom");
 for(let n=1;n<=6;n++){const D=mpDet(ord.slice(0,n).map(i=>ord.slice(0,n).map(j=>Q[i][j]))),stats=mpStats(D.n);minors.push({n,stats,fraction:mpRFSerialize(D)});if(!wrongCenter&&(stats.negative||!stats.strictEdge))throw Error("P6 coefficient positivity failure "+ax+","+gi+","+n);}
 return{ax,gi,wrongCenter,ms:Date.now()-mpStart,ops:mpOps,maxTerms:mpMaxTerms,entryStats:Q.map(row=>row.map(x=>({terms:x.n.size,degree:mpDegree(x.n)}))),data:Q.map(row=>row.map(mpRFSerialize)),atomData:mpAtoms.map(x=>({p:mpSerialize(x.p),label:x.label,positive:x.positive})),minors};
}

return mpP6Build(ax,gi);
}

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


function sfQP6(t,a,s,F,alpha,B,gamma,H=sfq(3,4),centerMultiplier=sfq(1)){
 const K=t.add(sfq(1,2)),k=K.sub(a.div(2)),delta=k.add(a),theta=a.div(t),env=t.sub(a),u=t.mul(2).sub(a),v=k.sub(sfq(1,2)),D0=sfq(2).add(a.div(delta)),C0=sfq(2).add(a.add(1).div(delta)),d=K.add(s),e=k.add(s),f=sfq(1).div(d),g=sfq(1).div(e),om=s.add(sfq(1,4)),VF=sfq(1).div(K.sq().mul(K.add(1).sq()).mul(4)),VG=sfq(1).div(k.sq().mul(k.add(1).sq()).mul(4)),sqrtVF=sfq(1).div(K.mul(K.add(1)).mul(2)),sqrtVG=sfq(1).div(k.mul(k.add(1)).mul(2));
 const P=Array.from({length:6},()=>Array.from({length:6},()=>sfq(0))),add=(i,j,z)=>{P[i][j]=P[i][j].add(z);if(i!==j)P[j][i]=P[j][i].add(z);},outer=(co,x)=>{for(let i=0;i<6;i++)for(let j=0;j<6;j++)P[i][j]=P[i][j].add(co.mul(x[i]).mul(x[j]));};
 P[0][0]=s.add(sfq(1,2)).add(alpha.mul(s).mul(2)).mul(2).div(d).sub(t.mul(om).mul(4).div(d.sq()));
 P[1][1]=s.add(sfq(1,2)).sub(alpha.mul(s).mul(2)).mul(2).div(d);
 P[2][2]=sfq(1).sub(a.mul(om).div(e.sq().mul(4))).mul(t).div(a);
 add(0,2,gamma.mul(C0).mul(s).add(s.add(sfq(1,2))).sub(t.mul(om).mul(2).div(d)).div(e.mul(2)));
 add(1,2,gamma.mul(C0).mul(s).sub(s.add(sfq(1,2))).div(e.mul(2)));
 P[3][3]=sfq(1).sub(t.mul(alpha).mul(f).mul(2)).sub(t.div(K.sq().mul(2)));
 add(3,4,sfq(1).sub(sfq(1).div(delta.mul(2))).add(gamma.mul(a).div(2).mul(D0.mul(g).add(sfq(1).div(delta.sq())))).sub(t.div(k)).div(a.mul(2)));
 P[4][4]=sfq(2).sub(K.sq().div(k.sq().mul(2))).mul(t).div(a.sq());
 add(4,5,u.div(a).mul(sfq(2).sub(K.sq().div(k.sq().mul(2)))).sub(K.div(k.sq().mul(4))).mul(t).div(a));
 P[5][5]=sfq(1).add(env.mul(t.mul(3).sub(a)).mul(2).div(a)).mul(t).div(a);
 add(3,5,t.div(K.mul(a).mul(2)).sub(env.div(a.mul(2)).mul(a.add(1).div(delta).sub(gamma.mul(a).mul(sfq(1).div(delta.sq()).add(D0.div(delta)))))));
 const gammaMax=sfq(1).div(D0.add(sfq(1).div(k.add(1)))),cmax=a.mul(gammaMax).div(2).mul(a.div(k.mul(delta)).add(sfq(1).div(k.mul(k.add(1))))),epsEll=cmax.mul(sqrtVG),cEll=sfq(1).sub(theta).div(4),cE=cEll.mul(sfq(1).add(epsEll)).mul(t).div(a);
 function source(z){
   const ff=sfq(1).div(K.add(z)),gg=sfq(1).div(k.add(z)),oo=z.add(sfq(1,4));
   const WA=[[t.mul(alpha).mul(ff).mul(-4),gamma.mul(t.mul(2).add(a.mul(2)).add(1)).mul(gg).div(delta.mul(2)),gamma.mul(env).mul(gg).div(delta.mul(-2))],[sfq(0),gamma.mul(gg).div(delta.mul(2)),gamma.mul(env).mul(gg).div(delta.mul(2))],[gamma.mul(C0).mul(v).mul(gg).div(-2),sfq(0),sfq(0)]];
   const WB=[[t.sq().mul(2).div(K).mul(ff).mul(oo.mul(ff).add(F.div(4))),t.div(k.mul(-2)).mul(ff).mul(oo.mul(gg).add(B.div(4))),sfq(0)],[sfq(0),sfq(0),sfq(0)],[t.sq().div(K.mul(2)).mul(gg).mul(oo.mul(ff).add(F.div(4))),t.div(k.mul(-8)).mul(gg).mul(oo.mul(gg).add(B.div(4))),sfq(0)]];
   const ae=sfq(1).add(a.mul(gg).div(2)),be=sfq(1).sub(a.mul(gamma).div(delta)).add(a.mul(gamma).mul(gg)),mc=[sfq(0),sfq(1).div(delta),t.mul(3).add(1).div(delta)];
   return WA.map((row,i)=>row.map((x,j)=>x.sub(WB[i][j]).sub(i===2?cE.mul(ae).mul(be).mul(mc[j]):sfq(0))));
 }
 const W=source(s),W0=source(sfq(0));for(let i=0;i<3;i++)for(let j=0;j<3;j++)add(i,3+j,W[i][j].sub(W0[i][j].mul(centerMultiplier)).add(i===0&&j===0?t.div(t.add(1).sq().mul(t.add(1)).mul(4)):sfq(0)));
 const bs=[sfq(0),sfq(0),sfq(0),t.div(K).mul(f),g.div(k.mul(-4)),sfq(0)],bbs=[sfq(0),sfq(0),sfq(0),t.div(K).mul(F),B.div(k.mul(-4)),sfq(0)];outer(t.mul(om).neg(),bs);outer(t.div(-4),bbs);
 const zLoss=a.mul(VG).div(16).add(a.div(t).mul(K.add(1).sq()).mul(VG).mul(gamma.mul(v).mul(sqrtVG).mul(2).add(t.mul(sqrtVF)).sq())).add(gamma.sq().mul(a).div(t).mul(v.sq()).mul(K.add(1)).mul(VG.sq()).mul(4));
 P[0][0]=P[0][0].sub(t.mul(sfq(1).sub(alpha.mul(2))).div(K.sq().mul(4)).add(sfq(1,16)).mul(f.sq()));P[1][1]=P[1][1].sub(f.div(16));P[2][2]=P[2][2].sub(zLoss.mul(t).div(a));
 outer(t.mul(alpha.sq()).mul(VF).mul(-64),[s.mul(f),sfq(0),sfq(0),t.mul(f).div(-2),sfq(0),sfq(0)]);
 outer(t.mul(alpha.sq()).mul(VF).mul(-32),[sfq(0),s.mul(f),sfq(0),sfq(0),sfq(0),sfq(0)]);
 outer(t.mul(gamma.sq()).mul(VG).mul(-4),[sfq(0),sfq(0),s.mul(g),sfq(0),g.div(2),sfq(0)]);
 const ae=sfq(1).add(a.mul(g).div(2)),be=sfq(1).sub(a.mul(gamma).div(delta)).add(a.mul(gamma).mul(g)),mc=[sfq(0),sfq(1).div(delta),t.mul(3).add(1).div(delta)];
 outer(cE.neg(),[sfq(0),sfq(0),ae,sfq(0),be.div(delta),be.mul(t.mul(3).add(1)).div(delta)]);
 for(let j=0;j<3;j++)add(2,3+j,cE.mul(ae).mul(be).mul(mc[j]));
 P[2][2]=P[2][2].sub(cEll.mul(epsEll.add(epsEll.sq())).mul(t).div(a));
 const co=H.mul(sfq(1).sub(theta));P[2][2]=P[2][2].sub(co.mul(t).div(a));
 outer(co.mul(t).div(a.sq()).neg(),[sfq(0),sfq(0),sfq(0),sfq(0),sfq(1),u]);P[5][5]=P[5][5].sub(co.mul(t.sq()).div(a.mul(t.add(1))));
 return{P,K,k,delta,epsEll};
}

const sfQZero=n=>Array.from({length:n},()=>Array.from({length:n},()=>sfq(0)));
const sfQDiag=v=>v.map((x,i)=>v.map((_,j)=>i===j?x:sfq(0)));
const sfQSub=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.sub(B[i][j])));
const sfQAdd=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.add(B[i][j])));
const sfQScale=(A,c)=>A.map(row=>row.map(x=>x.mul(c)));
const sfQCong=(G,A)=>sfQM(sfQT(G),sfQM(A,G));
function sfQP6FormMatrix(n,fn){const b=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>sfq(i===j?1:0))),v=b.map(fn),M=sfQDiag(v);for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)M[i][j]=M[j][i]=fn(b[i].map((z,k)=>z.add(b[j][k]))).sub(v[i]).sub(v[j]).div(2);return M;}
function sfQP6IntegratedGuard(weights,eps,x,mu,H=sfq(3,4)){
 const t=sfdot(weights,weights),a=weights[0].sq(),w=weights[0],K=t.add(sfq(1,2)),k=K.sub(a.div(2)),delta=k.add(a);
 if(!eps.sq().eq(t))throw Error("choose rational square total for original-energy guard");
 const F=sfdot(mu,x.map(s=>sfq(1).div(K.add(s)))),B=sfdot(mu,x.map(s=>sfq(1).div(k.add(s)))),alpha=sfq(1).div(sfq(2).add(F)),gamma=sfq(1).div(sfq(2).add(a.div(delta)).add(B)),P=x.map(s=>sfQP6(t,a,s,F,alpha,B,gamma,H).P),u=x.length;
 const form=(v,which)=>{
  const field=j=>{const f=v.slice(j*(u-1),(j+1)*(u-1));f.push(sfdot(f,mu.slice(0,u-1)).neg().div(mu[u-1]));return f;},L=field(0),D=field(1),ZZ=field(2),[pp,rr,hh]=v.slice(3*(u-1)),p=pp.div(eps),r=rr.mul(eps).div(w),h=hh.mul(eps).div(w),Z0=ZZ.map(z=>z.mul(eps).div(w)),A=L.map((z,i)=>z.add(D[i]).div(eps)),C=L.map((z,i)=>z.sub(D[i]).div(eps).add(p)),Z=Z0.map(z=>z.add(h)),rho=r.add(t.mul(2).sub(a).mul(h)).div(w);
  if(which==="lower")return sfdot(mu,x.map((_,i)=>{const v=[L[i],D[i],ZZ[i],pp,rr,hh];return sfdot(v,sfMatVec(P[i],v));}));
  const exact=sfQTriangleRemainder(weights,x,mu,A,C,[Z,x.map(()=>sfq(0))],[rho,sfq(0)]).remainder.sub(H.mul(sfq(1).sub(a.div(t))).mul(rho.sq().add(sfdot(mu,Z0.map(z=>z.sq()))).add(t.mul(h.sq()).div(t.add(1)))));
  return exact;
 };
 const E=sfQP6FormMatrix(3*u,v=>form(v,"actual")),lower=sfQP6FormMatrix(3*u,v=>form(v,"lower")),gap=sfQSub(E,lower),I=sfQInertia(gap),IL=sfQInertia(lower);if(I.negative)throw Error("P6 integrated lower bound failed");
 return{sigma:t.str(),a:a.str(),u,dimension:3*u,gap:I,lower:IL};
}
function sfQP6EqualM1Guard(t,s,F,alpha,B,gamma){
 // sigma=t=a; QP6 is in the root-cleared coordinates.
 const one=sfq(1),P6=sfQP6(t,t,s,F,alpha,B,gamma).P,wt=one;
 // Compare to the known m1 P5 with rational sqrt(t), supplied separately in wrapper.
 return P6;
}

function sfQP6Tri(pole,c){const lo=sfq(1).div(pole.add(1)),hi=sfq(1).div(pole),left=c.add(lo),right=c.add(hi);return[[lo,sfq(1).div(left)],[hi,sfq(1).div(right)],[left.mul(right).mul(2).div(left.add(right)).sub(c),sfq(2).div(left.add(right))]];}
function sfQP6Bary(V,x,y){const A=V[0],B=V[1],C=V[2],dx=x.sub(C[0]),dy=y.sub(C[1]),aa=A[0].sub(C[0]),bb=B[0].sub(C[0]),cc=A[1].sub(C[1]),dd=B[1].sub(C[1]),det=aa.mul(dd).sub(bb.mul(cc)),r=dx.mul(dd).sub(bb.mul(dy)).div(det),s=aa.mul(dy).sub(dx.mul(cc)).div(det);return[r,s,sfq(1).sub(r).sub(s)];}
function sfQP6ConcavityGuard(t,a,x,mu,s){
 const K=t.add(sfq(1,2)),k=K.sub(a.div(2)),delta=k.add(a),D0=sfq(2).add(a.div(delta)),F=sfdot(mu,x.map(v=>sfq(1).div(K.add(v)))),B=sfdot(mu,x.map(v=>sfq(1).div(k.add(v)))),alpha=sfq(1).div(sfq(2).add(F)),gamma=sfq(1).div(D0.add(B)),A=sfQP6Tri(K,sfq(2)),G=sfQP6Tri(k,D0),wa=sfQP6Bary(A,F,alpha),wg=sfQP6Bary(G,B,gamma),P=sfQP6(t,a,s,F,alpha,B,gamma).P;let mix=sfQZero(6);
 if([...wa,...wg].some(v=>!v.ge(0)))throw Error("P6 triangle enclosure");
 for(let ax=0;ax<3;ax++)for(let gi=0;gi<3;gi++)mix=sfQAdd(mix,sfQScale(sfQP6(t,a,s,A[ax][0],A[ax][1],G[gi][0],G[gi][1]).P,wa[ax].mul(wg[gi])));
 const I=sfQInertia(sfQSub(P,mix));if(I.negative)throw Error("P6 concavity");return{sigma:t.str(),a:a.str(),s:s.str(),gap:I};
}
function sfQP6WrongCenterGuard(){
 const t=sfq(10000),a=t,s=sfq(0),K=t.add(sfq(1,2)),k=t.add(1).div(2),D0=sfq(2).add(t.div(t.mul(3).add(1).div(2))),A=sfQP6Tri(K,sfq(2)),G=sfQP6Tri(k,D0),P=sfQP6(t,a,s,A[0][0],A[0][1],G[0][0],G[0][1],sfq(1,2),sfq(2)).P,I=sfQInertia(P);if(!I.negative)throw Error("P6 wrong center did not fail");return{sigma:t.str(),a:a.str(),s:s.str(),inertia:I};
}

function sfQP6ExactDecomposition(wt,eps,x,mu,L,D,Z0,p,r,h){
 const t=sfdot(wt,wt),a=wt[0].sq(),w=wt[0],K=t.add(sfq(1,2)),k=K.sub(a.div(2)),delta=k.add(a),D0=sfq(2).add(a.div(delta)),C0coef=sfq(2).add(a.add(1).div(delta)),v=k.sub(sfq(1,2)),env=t.sub(a),u=t.mul(2).sub(a),avg=z=>sfdot(mu,z),f=x.map(s=>sfq(1).div(K.add(s))),g=x.map(s=>sfq(1).div(k.add(s))),F=avg(f),B=avg(g),alpha=sfq(1).div(sfq(2).add(F)),gamma=sfq(1).div(D0.add(B)),A=L.map((z,i)=>z.add(D[i]).div(eps)),C=L.map((z,i)=>z.sub(D[i]).div(eps).add(p)),Z=Z0.map(z=>z.add(h)),rho=r.add(u.mul(h)).div(w),direct=sfQTriangleRemainder(wt,x,mu,A,C,[Z,x.map(()=>sfq(0))],[rho,sfq(0)]);
 let local=sfq(0);for(let i=0;i<x.length;i++){
 const s=x[i],d=K.add(s),e=k.add(s),ff=f[i],gg=g[i],om=s.add(sfq(1,4)),M0LL=s.add(sfq(1,2)).add(alpha.mul(s).mul(2)).mul(2).mul(ff).sub(t.mul(om).mul(ff.sq()).mul(4)),M0DD=s.add(sfq(1,2)).sub(alpha.mul(s).mul(2)).mul(2).mul(ff),M0ZZ=sfq(1).sub(a.mul(om).mul(gg.sq()).div(4)),M0LZ=w.div(eps.mul(2)).mul(gg).mul(gamma.mul(C0coef).mul(s).add(s.add(sfq(1,2))).sub(t.mul(om).mul(ff).mul(2))),M0DZ=w.div(eps.mul(2)).mul(gg).mul(gamma.mul(C0coef).mul(s).sub(s.add(sfq(1,2))));
 const WA=[[t.mul(eps).mul(alpha).mul(ff).mul(-4),w.mul(gamma).mul(t.mul(2).add(a.mul(2)).add(1)).mul(gg).div(eps.mul(delta).mul(2)),w.mul(gamma).mul(env).mul(gg).div(eps.mul(delta).mul(-2))],[sfq(0),w.mul(gamma).mul(gg).div(eps.mul(delta).mul(2)),w.mul(gamma).mul(env).mul(gg).div(eps.mul(delta).mul(2))],[w.mul(gamma).mul(C0coef).mul(v).mul(gg).div(-2),sfq(0),sfq(0)]];
 const af=[eps.mul(ff).mul(2),sfq(0),w.mul(gg).div(2)],bv=[t.sq().div(K).mul(ff),w.div(k.mul(-4)).mul(gg),sfq(0)],bb=[t.sq().div(K).mul(F),w.div(k.mul(-4)).mul(B),sfq(0)],ffields=[L[i],D[i],Z0[i]],means=[p,r,h],meanA=[[t.sub(t.sq().mul(alpha).mul(ff).mul(2)).sub(t.sq().div(K.sq().mul(2))),sfq(1).sub(sfq(1).div(delta.mul(2))).add(gamma.mul(a).div(2).mul(D0.mul(gg).add(sfq(1).div(delta.sq())))).sub(t.div(k)).div(w.mul(2)),t.div(K.mul(w).mul(2)).sub(env.div(w.mul(2)).mul(a.add(1).div(delta).sub(gamma.mul(a).mul(sfq(1).div(delta.sq()).add(D0.div(delta))))))],[sfq(0),sfq(2).sub(K.sq().div(k.sq().mul(2))).div(a),u.div(a).mul(sfq(2).sub(K.sq().div(k.sq().mul(2)))).sub(K.div(k.sq().mul(4)))],[sfq(0),sfq(0),sfq(1).add(env.mul(t.mul(3).sub(a)).mul(2).div(a))]];
 for(let j=0;j<3;j++)for(let l=j+1;l<3;l++)meanA[l][j]=meanA[j][l];
 let z=M0LL.mul(L[i].sq()).add(M0DD.mul(D[i].sq())).add(M0ZZ.mul(Z0[i].sq())).add(M0LZ.mul(L[i]).mul(Z0[i]).mul(2)).add(M0DZ.mul(D[i]).mul(Z0[i]).mul(2));
 for(let j=0;j<3;j++)for(let l=0;l<3;l++)z=z.add(ffields[j].mul(means[l]).mul(WA[j][l].sub(af[j].mul(om.mul(bv[l]).add(bb[l].div(4))))).mul(2));
 z=z.add(sfdot(means,sfMatVec(meanA,means))).sub(om.mul(sfdot(bv,means).sq())).sub(sfdot(bb,means).sq().div(4));local=local.add(mu[i].mul(z));
 }
 const Lf=avg(L.map((z,i)=>z.mul(f[i]))),Df=avg(D.map((z,i)=>z.mul(f[i]))),Lg=avg(L.map((z,i)=>z.mul(g[i]))),Dg=avg(D.map((z,i)=>z.mul(g[i]))),Zg=avg(Z0.map((z,i)=>z.mul(g[i]))),collect=t.mul(sfq(1).sub(alpha.mul(2))).mul(Lf.sq()).neg().add(t.mul(alpha).mul(Df.sq()).mul(2)).sub(eps.mul(w).mul(Lf).mul(Zg).div(2)).sub(a.mul(Zg.sq()).div(16)).add(gamma.mul(w).mul(v).div(eps).mul(Lg.sub(Dg)).mul(Zg)),M=r.add(t.mul(3).add(1).mul(h)).div(delta),ell=x.map((s,i)=>sfq(1).add(a.mul(g[i]).div(2)).mul(Z0[i]).neg().sub(sfq(1).sub(a.mul(gamma).div(delta)).add(a.mul(gamma).mul(g[i])).mul(M)).sub(a.mul(gamma).div(2).mul(sfq(1).div(delta).sub(g[i])).mul(Zg))),ellLoss=sfq(1).sub(a.div(t)).div(4).mul(avg(ell.map(z=>z.sq()))),reconstructed=local.add(collect).sub(ellLoss);
 if(!reconstructed.eq(direct.base.sub(direct.B)))throw Error("P6 exact A-B decomposition mismatch");
 return true;
}
function sfQP6DecompositionGuard(wt,eps,x,mu){
 const u=x.length,n=3*u,field=v=>{const z=[...v];z.push(sfdot(v,mu.slice(0,u-1)).neg().div(mu[u-1]));return z;},fn=v=>sfQP6ExactDecomposition(wt,eps,x,mu,field(v.slice(0,u-1)),field(v.slice(u-1,2*(u-1))),field(v.slice(2*(u-1),3*(u-1))),...v.slice(3*(u-1))),basis=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>sfq(i===j?1:0)));
 let count=0;for(let i=0;i<n;i++){fn(basis[i]);count++;for(let j=i+1;j<n;j++){fn(basis[i].map((z,k)=>z.add(basis[j][k])));count++;}}
 return{sigma:sfdot(wt,wt).str(),a:wt[0].sq().str(),u,coefficientIdentities:count};
}

function sfQMPPoly(poly,a,b,s){let z=sfq(0);const pow=(x,n)=>{let y=sfq(1);while(n){if(n&1)y=y.mul(x);n>>=1;if(n)x=x.sq();}return y;};const deg=[0,0,0];for(const[k]of poly){const ex=[Math.floor(k/65536),Math.floor(k/256)%256,k%256];for(let j=0;j<3;j++)deg[j]=Math.max(deg[j],ex[j]);}const pows=[a,b,s].map((v,i)=>{const p=[sfq(1)];for(let j=0;j<deg[i];j++)p.push(p.at(-1).mul(v));return p;});for(const[k,v]of poly)z=z.add(new SFQ(BigInt(v)).mul(pows[0][Math.floor(k/65536)]).mul(pows[1][Math.floor(k/256)%256]).mul(pows[2][k%256]));return z;}
function sfQMPFraction(x,atoms,a,b,s){let num=sfQMPPoly(x.n,a,b,s),den=new SFQ(BigInt(x.d));for(const[id,ex]of x.f){const v=atoms[id];for(let k=0;k<ex;k++)den=den.mul(v);}return num.div(den);}
function sfQP6CertificateGuard(data,points){
 let entries=0,determinants=0;
 for(const [a,b,s]of points){const t=a.add(b),K=t.add(sfq(1,2)),k=K.sub(a.div(2)),delta=k.add(a),tri=(pole,c)=>{const lo=sfq(1).div(pole.add(1)),hi=sfq(1).div(pole),left=c.add(lo),right=c.add(hi);return[[lo,sfq(1).div(left)],[hi,sfq(1).div(right)],[left.mul(right).mul(2).div(left.add(right)).sub(c),sfq(2).div(left.add(right))]];},A=tri(K,sfq(2)),G=tri(k,sfq(2).add(a.div(delta))),raw=sfQP6(t,a,s,A[data.ax][0],A[data.ax][1],G[data.gi][0],G[data.gi][1]).P,scale=[K.add(s),K.add(s),k.add(s),K.add(s),k.add(s),k.add(s)],direct=raw.map((row,i)=>row.map((z,j)=>z.mul(scale[i]).mul(scale[j]))),atomVals=data.atomData.map(x=>sfQMPPoly(x.p,a,b,s)),saved=data.data.map(row=>row.map(x=>sfQMPFraction(x,atomVals,a,b,s)));
 for(let i=0;i<6;i++)for(let j=0;j<6;j++){if(!saved[i][j].eq(direct[i][j]))throw Error("independent P6 polynomial entry mismatch");entries++;}
 const ord=[5,1,2,4,0,3];for(let n=1;n<=6;n++){const D=sfdet(ord.slice(0,n).map(i=>ord.slice(0,n).map(j=>direct[i][j]))),cert=sfQMPFraction(data.minors[n-1].fraction,atomVals,a,b,s);if(!D.eq(cert)||!D.gt(0))throw Error("independent P6 determinant mismatch");determinants++;}
 }
 return{ax:data.ax,gi:data.gi,points:points.length,entries,determinants,status:"PASS"};
}

function sfQMeanTwoSources(wt,x,mu,B){
 const t=wt.sq(),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),{F,F2,B:Bg,B2,a:alpha,gam,cov,vf,vg}=B.loc;
 return x.map(s=>{const f=sfq(1).div(K.add(s)),g=sfq(1).div(k.add(s)),w=s.add(sfq(1,4)),A=[[t.mul(wt).mul(alpha).mul(f).mul(-4),wt.mul(gam).mul(t.mul(4).add(1)).mul(g).div(delta.mul(2))],[sfq(0),wt.mul(gam).mul(g).div(delta.mul(2))],[wt.mul(t).mul(K).mul(gam).mul(g).neg().div(delta),sfq(0)]],Ba=[[wt.mul(t.sq()).mul(f).mul(w.mul(f).add(F.div(4))).mul(2).div(K),-1],[sfq(0),sfq(0)],[wt.mul(t.sq()).mul(g).mul(w.mul(f).add(F.div(4))).div(K.mul(2)),-1]],R=[[wt.mul(t.sq()).mul(alpha.sq()).mul(vf).mul(s).mul(f).mul(f.sub(F)).mul(-2),wt.mul(t).mul(alpha).mul(gam).mul(cov).mul(s).mul(f).mul(g.sub(Bg)).div(2)],[sfq(0),sfq(0)],[wt.mul(t.sq()).mul(alpha).mul(gam).mul(cov).mul(s).mul(g).mul(f.sub(F)).div(-2),wt.mul(t).mul(gam.sq()).mul(vg).mul(s).mul(g).mul(g.sub(Bg)).div(8)]];
 Ba[0][1]=wt.mul(t).mul(f).mul(w.mul(g).add(Bg.div(4))).neg().div(k.mul(2));Ba[2][1]=wt.mul(t).mul(g).mul(w.mul(g).add(Bg.div(4))).neg().div(k.mul(8));
 return sfQSub(sfQSub(A,Ba),sfQScale(R,sfq(8)));
 });
}

function sfQUniversalMeanPoint(wt,s,F,alpha,B,gam){
 const t=wt.sq(),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),d=K.add(s),e=k.add(s),f=sfq(1).div(d),g=sfq(1).div(e),w=s.add(sfq(1,4)),vf=sfq(1).div(K.sq().mul(K.add(1).sq()).mul(4)),vg=sfq(1).div(k.sq().mul(k.add(1).sq()).mul(4)),ra=sfq(1).sub(alpha.mul(2)),u=K.mul(gam).mul(2).div(delta).sub(sfq(1,2)),bL=K.add(1).mul(gam).div(k.mul(k.add(1)).mul(2)).add(sfq(1).div(K.mul(2))).sq(),bD=K.add(1).mul(vg),zLoss=t.mul(vg).div(16).add(t.sq().mul(vg).mul(bL.add(gam.sq().mul(bD)))),LL=s.add(sfq(1,2)).add(alpha.mul(s).mul(2)).mul(2).div(d).sub(t.mul(w).mul(4).div(d.sq())).sub(t.mul(ra).div(K.sq().mul(4)).add(sfq(1,16)).div(d.sq())),DD=s.add(sfq(1,2)).sub(alpha.mul(s).mul(2)).mul(2).div(d).sub(sfq(1).div(d.mul(16))),ZZ=sfq(1).sub(t.mul(w).div(e.sq().mul(4))).sub(zLoss),cp=u.mul(s).add(w.mul(s.add(sfq(1,2))).div(d)).div(e),cm=u.mul(s).sub(sfq(1,4)).div(e),P=sfQZero(5);
 P[0][0]=LL;P[1][1]=DD;P[2][2]=ZZ;P[0][2]=P[2][0]=cp;P[1][2]=P[2][1]=cm;
 const loc={F,F2:F.sq(),B,B2:B.sq(),a:alpha,gam,vf:sfq(0),vg:sfq(0),cov:sfq(0)},sources=sfQMeanTwoSources(wt,[sfq(0),s],[],{loc}),cor=wt.mul(t).div(t.add(1).mul(t.add(1)).mul(t.add(1)).mul(4));
 for(let a=0;a<3;a++)for(let b=0;b<2;b++){let v=sources[1][a][b].sub(sources[0][a][b]);if(a===0&&b===0)v=v.add(cor);P[a][3+b]=P[3+b][a]=v;}
 const D0=sfq(2).add(t.div(delta)),a12=sfq(1).sub(sfq(1).div(delta.mul(2))).add(gam.mul(t).div(2).mul(D0.mul(g).add(sfq(1).div(delta.sq())))).sub(t.mul(2).div(t.add(1))).sub(t.div(K.mul(t.add(1)))).div(2);
 P[3][3]=t.sq().mul(sfq(1).div(K).sub(alpha.mul(f).mul(2))).sub(t.sq().sq().div(K.sq()).mul(w.mul(f.sq()).add(F.sq().div(4))));
 P[4][4]=K.add(1).div(t.add(1).sq()).sub(t.sq().div(k.sq().mul(16)).mul(w.mul(g.sq()).add(B.sq().div(4))));
 P[3][4]=P[4][3]=a12.add(t.sq().mul(t).div(K.mul(k).mul(4)).mul(w.mul(f).mul(g).add(F.mul(B).div(4))));
 const ranks=[
 {c:t.mul(alpha.sq()).mul(vf).mul(64),r:[s.mul(f),sfq(0),sfq(0),t.mul(wt).mul(f).div(-2),sfq(0)]},
 {c:t.mul(alpha.sq()).mul(vf).mul(32),r:[sfq(0),s.mul(f),sfq(0),sfq(0),sfq(0)]},
 {c:t.mul(gam.sq()).mul(vg).mul(4),r:[sfq(0),sfq(0),s.mul(g),sfq(0),wt.mul(g).div(2)]}
 ];
 for(const z of ranks)for(let a=0;a<5;a++)for(let b=0;b<5;b++)P[a][b]=P[a][b].sub(z.c.mul(z.r[a]).mul(z.r[b]));
 return{P,scaled:sfQCong(sfQDiag([d,d,e,d.div(wt),e.div(wt)]),P)};
}
function sfQTriangleInverse(lo,hi,shift){const a=shift.add(lo),b=shift.add(hi);return[{x:lo,y:sfq(1).div(a)},{x:hi,y:sfq(1).div(b)},{x:a.mul(b).mul(2).div(a.add(b)).sub(shift),y:sfq(2).div(a.add(b))}];}
function sfQEvalT(poly,t){let v=sfq(0);for(let j=poly.length-1;j>=0;j--)v=v.mul(t).add(new SFQ(BigInt(poly[j])));return v;}
function sfQEvalTS(poly,t,s){let v=sfq(0);for(let j=poly.length-1;j>=0;j--)v=v.mul(s).add(sfQEvalT(poly[j],t));return v;}

function sfP6AuxiliaryChecks(){
 const formalX=[sfq(0),sfq(1,2),sfq(1)],formalMu=[sfq(49,100),sfq(1,50),sfq(49,100)],physicalX=[sfq(1,10),sfq(1,5),sfq(7,10)],profiles=[
 [[sfq(3,5),sfq(4,5)],sfq(1),formalX,formalMu],
 [[sfq(4,5),sfq(3,5)],sfq(1),physicalX,physicalX],
 [[sfq(3,2),sfq(2)],sfq(5,2),[sfq(0),sfq(1,5),sfq(7,10),sfq(1)],[sfq(1,10),sfq(1,5),sfq(3,10),sfq(2,5)]],
 [[sfq(3,50),sfq(2,25)],sfq(1,10),formalX,formalMu],
 [[sfq(6),sfq(8)],sfq(10),formalX,formalMu],
 [[sfq(1),sfq(0)],sfq(1),physicalX,physicalX]
 ],integrated=profiles.map(p=>sfQP6IntegratedGuard(...p)),decomposition=profiles.map(p=>sfQP6DecompositionGuard(...p)),concavity=[];
 for(const [weights,eps,x,mu]of profiles)for(const s of [sfq(0),sfq(1,3),sfq(1)])concavity.push(sfQP6ConcavityGuard(eps.sq(),weights[0].sq(),x,mu,s));
 let m1Identities=0;
 for(const wt of [sfq(1,10),sfq(1,2),sfq(1),sfq(3),sfq(10)])for(const s of [sfq(0),sfq(1,3),sfq(1)]){
  const t=wt.sq(),K=t.add(sfq(1,2)),k=t.add(1).div(2),D0=sfq(2).add(t.div(t.mul(3).add(1).div(2))),TA=sfQTriangleInverse(sfq(1).div(K.add(1)),sfq(1).div(K),sfq(2)),TG=sfQTriangleInverse(sfq(1).div(k.add(1)),sfq(1).div(k),D0);
  for(let ax=0;ax<3;ax++)for(let gi=0;gi<3;gi++){
   const P5=sfQUniversalMeanPoint(wt,s,TA[ax].x,TA[ax].y,TG[gi].x,TG[gi].y).P,G=Array.from({length:5},()=>Array(6).fill(sfq(0)));for(let i=0;i<3;i++)G[i][i]=sfq(1);G[3][3]=sfq(1).div(wt);G[4][4]=sfq(1).div(wt);
   const expected=sfQCong(G,P5),co=[sfq(0),sfq(0),sfq(0),sfq(1).div(K.mul(2)),sfq(1).div(t.add(1)),sfq(1)],P6=sfQP6(t,t,s,TA[ax].x,TA[ax].y,TG[gi].x,TG[gi].y).P;
   for(let i=0;i<6;i++)for(let j=0;j<6;j++){if(!P6[i][j].eq(expected[i][j].add(co[i].mul(co[j]))))throw Error("P6 m1 reduction");m1Identities++;}
  }
 }
 return{status:"PASS",integrated,decomposition,concavity,m1Identities,wrongCenter:sfQP6WrongCenterGuard()};
}
function sfP6Verify(options={}){
 const families=options.families||Array.from({length:9},(_,j)=>[Math.floor(j/3),j%3]),summaries=[],guards=[];
 for(const [ax,gi]of families){
  const d=sfP6BuildCertificate(ax,gi),points=[[sfq(1,100),sfq(0),sfq(0)],[sfq(9,25),sfq(16,25),sfq(1,3)],[sfq(7,3),sfq(2,7),sfq(1)],[sfq(100),sfq(1000),sfq(2)]];
  guards.push(sfQP6CertificateGuard(d,points));summaries.push({ax,gi,minors:d.minors.map(m=>({n:m.n,...m.stats})),ms:d.ms});
 }
 return{status:"PASS",families:summaries,leadingMinors:families.length*6,positiveCoefficients:summaries.reduce((v,x)=>v+x.minors.reduce((u,m)=>u+m.positive,0),0),entryGuards:guards.reduce((v,x)=>v+x.entries,0),determinantGuards:guards.reduce((v,x)=>v+x.determinants,0),auxiliary:options.skipAuxiliary?null:sfP6AuxiliaryChecks()};
}

