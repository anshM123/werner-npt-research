// Certificate of the energy theorem for one copy carrying the letter 2 (NOTES.md, Sections 2 and 3).
// A function body that returns its verdict: run it with `python run_js.py one_letter_certificate.js`.
const CERTS=(()=>{

let rfOps=0,rfMaxDegree=0,rfStart=Date.now(),rfDeadline=rfStart+45000;
function riGcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b){const t=a%b;a=b;b=t;}return a;}
function rpTrim(a){while(a.length>1&&a.at(-1)===0n)a.pop();return a;}
const rpZero=a=>a.length===1&&a[0]===0n;
function rpContent(a){let d=0n;for(const x of a)d=riGcd(d,x);return d;}
function rpPrimitive(a){a=rpTrim(a.slice());if(rpZero(a))return[0n];let c=rpContent(a);if(a.at(-1)<0n)c=-c;return a.map(x=>x/c);}
function rpAdd(a,b){return rpTrim(Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]??0n)+(b[i]??0n)));}
const rpNeg=a=>a.map(x=>-x),rpSub=(a,b)=>rpAdd(a,rpNeg(b));
function rpMul(a,b){if(rpZero(a)||rpZero(b))return[0n];const c=Array(a.length+b.length-1).fill(0n);for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]+=a[i]*b[j];return rpTrim(c);}
function rpPseudoRem(a,b){let r=a.slice();while(!rpZero(r)&&r.length>=b.length){const h=r.length-b.length,lc=r.at(-1),bc=b.at(-1);r=r.map(x=>x*bc);for(let j=0;j<b.length;j++)r[h+j]-=lc*b[j];r=rpPrimitive(r);if(Date.now()>rfDeadline)throw Error("RF computation time cap");}return r;}
function rpGcd(a,b){
 if(rpZero(a))return b.slice();if(rpZero(b))return a.slice();
 const content=riGcd(rpContent(a),rpContent(b));let x=rpPrimitive(a),y=rpPrimitive(b);if(x.length<y.length)[x,y]=[y,x];if(y.length===1)return[content];
 while(!rpZero(y)){const z=rpPseudoRem(x,y);x=y;y=z;if(y.length===1&&!rpZero(y))return[content];}
 return rpPrimitive(x).map(v=>v*content);
}
function rpDiv(a,b){if(b.length===1){if(a.some(x=>x%b[0]!==0n))throw Error("nonexact scalar poly division");return a.map(x=>x/b[0]);}
 let r=a.slice(),q=Array(Math.max(1,a.length-b.length+1)).fill(0n);while(!rpZero(r)&&r.length>=b.length){const h=r.length-b.length,lc=r.at(-1);if(lc%b.at(-1)!==0n)throw Error("nonexact polynomial quotient");const z=lc/b.at(-1);q[h]=z;for(let j=0;j<b.length;j++)r[h+j]-=z*b[j];rpTrim(r);}if(!rpZero(r))throw Error("nonzero polynomial division remainder");return rpTrim(q);
}
function rfTick(a,b){rfOps++;rfMaxDegree=Math.max(rfMaxDegree,a.length-1,b.length-1);if(rfMaxDegree>150)throw Error("RF degree resource guard");if(Date.now()>rfDeadline)throw Error("RF computation time cap");}
class SFRF{
 constructor(n,d=[1n],reduced=false){n=rpTrim(n);d=rpTrim(d);if(rpZero(d))throw Error("zero RF denominator");if(rpZero(n)){this.n=[0n];this.d=[1n];return;}if(!reduced){const g=rpGcd(n,d);n=rpDiv(n,g);d=rpDiv(d,g);}if(d.at(-1)<0n){n=rpNeg(n);d=rpNeg(d);}this.n=n;this.d=d;rfTick(n,d);}
 static constant(a,b=1){return new SFRF([BigInt(a)],[BigInt(b)]);}
 add(y){y=sfq(y);if(rpZero(this.n))return y;if(rpZero(y.n))return this;const same=this.d.length===y.d.length&&this.d.every((v,i)=>v===y.d[i]);if(same)return new SFRF(rpAdd(this.n,y.n),this.d);const g=rpGcd(this.d,y.d),x=rpDiv(this.d,g),z=rpDiv(y.d,g);return new SFRF(rpAdd(rpMul(this.n,z),rpMul(y.n,x)),rpMul(this.d,z));}
 sub(y){return this.add(sfq(y).neg());}
 neg(){return new SFRF(rpNeg(this.n),this.d,true);}
 mul(y){y=sfq(y);if(rpZero(this.n)||rpZero(y.n))return sfq(0);const g1=rpGcd(this.n,y.d),g2=rpGcd(y.n,this.d);return new SFRF(rpMul(rpDiv(this.n,g1),rpDiv(y.n,g2)),rpMul(rpDiv(this.d,g2),rpDiv(y.d,g1)),true);}
 div(y){y=sfq(y);if(rpZero(y.n))throw Error("RF divide by zero");return this.mul(new SFRF(y.d,y.n,true));}
 sq(){return this.mul(this);}
 eq(y){y=sfq(y);return this.n.length===y.n.length&&this.d.length===y.d.length&&this.n.every((v,i)=>v===y.n[i])&&this.d.every((v,i)=>v===y.d[i]);}
 str(){return"RF["+this.n.join(",")+"]/["+this.d.join(",")+"]";}
 at(x){const peval=p=>{let y=0n;for(let i=p.length-1;i>=0;i--)y=y*BigInt(x)+p[i];return y;};return{n:peval(this.n),d:peval(this.d)};}
}
const sfq=(a,b=1)=>a instanceof SFRF?a:SFRF.constant(a,b),sfsum=a=>a.reduce((v,x)=>v.add(x),sfq(0)),sfdot=(a,b)=>sfsum(a.map((v,i)=>v.mul(b[i])));

let sfRFProfileMoments=null;
let rfKnownFactors=[],rfFactorCache=new Map(),rfGcdCache=new Map();
const rpKey=a=>a.join(",");
function rpTryDiv(a,b){if(rpZero(a))return[0n];if(a.length<b.length)return null;let r=a.slice(),q=Array(a.length-b.length+1).fill(0n);while(!rpZero(r)&&r.length>=b.length){const h=r.length-b.length,z=r.at(-1);if(z%b.at(-1)!==0n)return null;const v=z/b.at(-1);q[h]=v;for(let j=0;j<b.length;j++)r[h+j]-=v*b[j];rpTrim(r);}return rpZero(r)?rpTrim(q):null;}
function rpFactorKnown(a){
 const key=rpKey(a);if(rfFactorCache.has(key))return rfFactorCache.get(key);
 const content=rpContent(a);if(!content)return null;let r=rpPrimitive(a),exp=rfKnownFactors.map(()=>0);
 for(let i=0;i<rfKnownFactors.length;i++){let q;while(r.length>1&&(q=rpTryDiv(r,rfKnownFactors[i]))!==null){r=q;exp[i]++;}}
 const result=r.length===1?{content,exp}:null;if(rfFactorCache.size<5000)rfFactorCache.set(key,result);return result;
}
function modp(a,p){const r=a%p;return r<0?r+p:r;}
function mpTrim(a){while(a.length>1&&a.at(-1)===0)a.pop();return a;}
function mpInv(a,p){let x=1,y=0,b=p;while(b){const q=Math.floor(a/b);[a,b]=[b,a-q*b];[x,y]=[y,x-q*y];}return modp(x,p);}
function mpGcdDegree(a,b,p=65521){
 const A=a.map(v=>Number((v%BigInt(p)+BigInt(p))%BigInt(p))),B=b.map(v=>Number((v%BigInt(p)+BigInt(p))%BigInt(p)));
 if(A.at(-1)===0||B.at(-1)===0)return -1;let x=mpTrim(A),y=mpTrim(B);
 while(!(y.length===1&&y[0]===0)){if(y.length===1)return 0;const r=x.slice(),inv=mpInv(y.at(-1),p);while(!(r.length===1&&r[0]===0)&&r.length>=y.length){const h=r.length-y.length,c=modp(r.at(-1)*inv,p);for(let j=0;j<y.length;j++)r[h+j]=modp(r[h+j]-c*y[j],p);mpTrim(r);}x=y;y=r;}return x.length-1;
}
function rpGcdSlowPrimitive(a,b){let x=rpPrimitive(a),y=rpPrimitive(b);if(x.length<y.length)[x,y]=[y,x];while(!rpZero(y)){if(y.length===1)return[1n];const z=rpPseudoRem(x,y);x=y;y=z;}return rpPrimitive(x);}
function rpGcdFast(a,b){
 if(rpZero(a))return b.slice();if(rpZero(b))return a.slice();const content=riGcd(rpContent(a),rpContent(b));
 if(a.length===1||b.length===1)return[content];
 let aa=rpPrimitive(a),bb=rpPrimitive(b),keyA=rpKey(aa),keyB=rpKey(bb),key=keyA<keyB?keyA+"|"+keyB:keyB+"|"+keyA;
 if(rfGcdCache.has(key))return rfGcdCache.get(key).map(v=>v*content);
 const fa=rpFactorKnown(aa),fb=rpFactorKnown(bb);let gp=[1n];
 if(fa&&fb){for(let i=0;i<rfKnownFactors.length;i++)for(let j=0;j<Math.min(fa.exp[i],fb.exp[i]);j++)gp=rpMul(gp,rfKnownFactors[i]);}
 else{
 if(fa||fb){if(fa&&!fb)[aa,bb]=[bb,aa];const fact=fa||fb;for(let i=0;i<rfKnownFactors.length;i++)for(let j=0;j<fact.exp[i];j++){const q=rpTryDiv(aa,rfKnownFactors[i]);if(q===null)break;aa=q;bb=rpDiv(bb,rfKnownFactors[i]);gp=rpMul(gp,rfKnownFactors[i]);}}
 if(aa.length>1&&bb.length>1&&mpGcdDegree(aa,bb)!==0)gp=rpMul(gp,rpGcdSlowPrimitive(aa,bb));
 }
 if(rfGcdCache.size<5000)rfGcdCache.set(key,gp);return gp.map(v=>v*content);
}
rpGcd=rpGcdFast;
function sfRFRegisterBesselFactors(data){
 const t=new SFRF([0n,1n]),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),mom=sfRFProfileMoments??[sfq(1),sfq(49,100),sfq(347,1000),sfq(2641,10000),sfq(20723,100000)],d=[K,sfq(1)],e=[k,sfq(1)],avg=p=>sfsum(p.map((v,j)=>v.mul(mom[j]))),Z=avg(spProd(spPow(d,2),spPow(e,2))),Nf=avg(spProd(d,spPow(e,2))),Nb=avg(spProd(spPow(d,2),e)),G=Z.mul(2).add(Nf),L=delta.mul(2).add(t).mul(Z).add(delta.mul(Nb));
 rfKnownFactors=[[0n,1n],[1n,1n],[1n,2n],[1n,3n],rpPrimitive(Z.n),rpPrimitive(G.n),rpPrimitive(L.n)];rfFactorCache.clear();rfGcdCache.clear();
 for(let i=0;i<rfKnownFactors.length;i++)for(let j=i+1;j<rfKnownFactors.length;j++)if(mpGcdDegree(rfKnownFactors[i],rfKnownFactors[j])!==0)throw Error("registered factors not certified pairwise coprime");
 return rfKnownFactors.map(f=>f.map(String));
}

const psC=c=>[sfq(c)],psZ=()=>psC(0),psS=()=>[sfq(0),sfq(1)];
const psTrim=a=>{while(a.length>1&&rpZero(a.at(-1).n))a.pop();return a;};
const psAdd=(a,b)=>psTrim(Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]??sfq(0)).add(b[i]??sfq(0))));
const psNeg=a=>a.map(v=>v.neg()),psSub=(a,b)=>psAdd(a,psNeg(b)),psScale=(a,c)=>a.map(v=>v.mul(c));
const psMul=(a,b)=>{const c=Array.from({length:a.length+b.length-1},()=>sfq(0));for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]=c[i+j].add(a[i].mul(b[j]));return psTrim(c);};
const psSq=a=>psMul(a,a),psSum=(...a)=>a.reduce(psAdd,psZ());
function sfRFUniversalMeanPoint(fi,gi,corFactor=sfq(1)){
 rfDeadline=Date.now()+45000;rfMaxDegree=0;rfOps=0;rfKnownFactors=[];rfGcdCache.clear();rfFactorCache.clear();
 const t=new SFRF([0n,1n]),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),triangle=(lo,hi,shift)=>{const a=shift.add(lo),b=shift.add(hi);return[{x:lo,y:sfq(1).div(a)},{x:hi,y:sfq(1).div(b)},{x:a.mul(b).mul(2).div(a.add(b)).sub(shift),y:sfq(2).div(a.add(b))}];},
 A=triangle(sfq(1).div(K.add(1)),sfq(1).div(K),sfq(2))[fi],G=triangle(sfq(1).div(k.add(1)),sfq(1).div(k),sfq(2).add(t.div(delta)))[gi],F=A.x,alpha=A.y,B=G.x,gam=G.y,D=[K,sfq(1)],E=[k,sfq(1)],S=psS(),S2=psSq(S),W=[sfq(1,4),sfq(1)],H=[sfq(1,2),sfq(1)],D2=psSq(D),E2=psSq(E),DE=psMul(D,E),vf=sfq(1).div(K.sq().mul(K.add(1).sq()).mul(4)),vg=sfq(1).div(k.sq().mul(k.add(1).sq()).mul(4)),rAlpha=sfq(1).sub(alpha.mul(2)),u=K.mul(gam).mul(2).div(delta).sub(sfq(1,2)),bL=K.add(1).mul(gam).div(k.mul(k.add(1)).mul(2)).add(sfq(1).div(K.mul(2))).sq(),bD=K.add(1).mul(vg),zLoss=t.mul(vg).div(16).add(t.sq().mul(vg).mul(bL.add(gam.sq().mul(bD)))),Q=Array.from({length:5},()=>Array.from({length:5},psZ));
 Q[0][0]=psSum(psScale(psMul(psAdd(H,psScale(S,alpha.mul(2))),D),sfq(2)),psScale(W,t.mul(-4)),psC(t.mul(rAlpha).div(K.sq().mul(4)).add(sfq(1,16)).neg()),psScale(S2,t.mul(alpha.sq()).mul(vf).mul(-64)));
 Q[1][1]=psSum(psScale(psMul(psSub(H,psScale(S,alpha.mul(2))),D),sfq(2)),psScale(D,sfq(-1,16)),psScale(S2,t.mul(alpha.sq()).mul(vf).mul(-32)));
 Q[2][2]=psSum(psScale(E2,sfq(1).sub(zLoss)),psScale(W,t.div(-4)),psScale(S2,t.mul(gam.sq()).mul(vg).mul(-4)));
 Q[0][2]=psAdd(psScale(psMul(S,D),u),psMul(W,H));Q[1][2]=psMul(psAdd(psScale(S,u),psC(sfq(-1,4))),D);
 Q[3][3]=psSum(psScale(D2,t.div(K)),psScale(D,t.mul(alpha).mul(-2)),psScale(psAdd(W,psScale(D2,F.sq().div(4))),t.sq().mul(t).neg().div(K.sq())),psC(t.sq().mul(t).mul(alpha.sq()).mul(vf).mul(-16)));
 Q[4][4]=psSum(psScale(E2,K.add(1).div(t.mul(t.add(1).sq()))),psScale(psAdd(W,psScale(E2,B.sq().div(4))),t.neg().div(k.sq().mul(16))),psC(t.mul(gam.sq()).mul(vg).neg()));
 const ac=sfq(1).sub(sfq(1).div(delta.mul(2))).add(gam.mul(t).div(delta.sq().mul(2))).sub(t.mul(2).div(t.add(1))).sub(t.div(K.mul(t.add(1)))).div(2),D0=sfq(2).add(t.div(delta));
 Q[3][4]=psSum(psScale(DE,ac.div(t)),psScale(D,gam.mul(D0).div(4)),psScale(psAdd(W,psScale(DE,F.mul(B).div(4))),t.sq().div(K.mul(k).mul(4))));
 const lp0=t.mul(alpha).mul(-4).div(K).sub(t.sq().div(K.sq().mul(K).mul(2))).sub(t.sq().mul(F).div(K.sq().mul(2))),lq0=gam.mul(t.mul(4).add(1)).div(delta.mul(k).mul(2)).add(t.div(k.mul(K).mul(2)).mul(sfq(1).div(k.mul(4)).add(B.div(4)))),zp0=t.mul(K).mul(gam).neg().div(delta.mul(k)).sub(t.sq().div(K.mul(k).mul(2)).mul(sfq(1).div(K.mul(4)).add(F.div(4)))),zq0=t.div(k.sq().mul(8)).mul(sfq(1).div(k.mul(4)).add(B.div(4))),cor=t.div(t.add(1).mul(t.add(1)).mul(t.add(1)).mul(4)).mul(corFactor);
 Q[0][3]=psSum(psScale(D,t.mul(alpha).mul(-4)),psScale(psAdd(W,psScale(D,F.div(4))),t.sq().mul(-2).div(K)),psScale(D2,cor.sub(lp0)),psScale(S,t.sq().mul(alpha.sq()).mul(vf).mul(32)));
 Q[0][4]=psSum(psScale(D,gam.mul(t.mul(4).add(1)).div(delta.mul(2))),psScale(psAdd(W,psScale(E,B.div(4))),t.div(k.mul(2))),psScale(DE,lq0.neg()));
 Q[1][4]=psScale(psMul(D,S),gam.neg().div(delta.mul(k).mul(2)));
 Q[2][3]=psSum(psScale(D,t.mul(K).mul(gam).neg().div(delta)),psScale(psAdd(W,psScale(D,F.div(4))),t.sq().neg().div(K.mul(2))),psScale(DE,zp0.neg()));
 Q[2][4]=psSum(psScale(psAdd(W,psScale(E,B.div(4))),t.div(k.mul(8))),psScale(E2,zq0.neg()),psScale(S,t.mul(gam.sq()).mul(vg).mul(-2)));
 for(let i=0;i<5;i++)for(let j=i+1;j<5;j++)Q[j][i]=Q[i][j];
 let Den=[1n];for(const row of Q)for(const p of row)for(const c of p)Den=rpMul(rpDiv(Den,rpGcd(Den,c.d)),c.d);
 const P=Q.map(row=>row.map(p=>p.map(c=>rpMul(c.n,rpDiv(Den,c.d)))));
 if(Den.some(v=>v<0n)||rpZero(Den))throw Error("denominator positivity certificate");
 return{fi,gi,P,Den,parameters:{F:{n:F.n,d:F.d},alpha:{n:alpha.n,d:alpha.d},B:{n:B.n,d:B.d},gam:{n:gam.n,d:gam.d}},maxDegreeT:Math.max(...P.flatMap(r=>r.flatMap(p=>p.map(c=>c.length-1)))),denDegree:Den.length-1,ops:rfOps};
}

const biTrim=A=>{while(A.length>1&&rpZero(A.at(-1)))A.pop();return A;};
const biZero=A=>A.length===1&&rpZero(A[0]);
const biAdd=(A,B)=>biTrim(Array.from({length:Math.max(A.length,B.length)},(_,i)=>rpAdd(A[i]??[0n],B[i]??[0n])));
const biNeg=A=>A.map(rpNeg),biSub=(A,B)=>biAdd(A,biNeg(B));
function biMul(A,B){if(biZero(A)||biZero(B))return[[0n]];const C=Array.from({length:A.length+B.length-1},()=>[0n]);for(let i=0;i<A.length;i++)for(let j=0;j<B.length;j++)C[i+j]=rpAdd(C[i+j],rpMul(A[i],B[j]));return biTrim(C);}
function biDet(A){const n=A.length,dp=Array(1<<n);dp[0]=[[1n]];for(let mask=0;mask<(1<<n);mask++)if(dp[mask]){let row=0;for(let j=0;j<n;j++)if(mask&(1<<j))row++;if(row===n)continue;for(let j=0;j<n;j++)if(!(mask&(1<<j))){let parity=0;for(let k=j+1;k<n;k++)if(mask&(1<<k))parity++;let term=biMul(dp[mask],A[row][j]);if(parity%2)term=biNeg(term);const m=mask|(1<<j);dp[m]=dp[m]?biAdd(dp[m],term):term;}}return dp[(1<<n)-1];}
function biStats(A){let pos=0,neg=0,zero=0;for(const r of A)for(const v of r)if(v>0n)pos++;else if(v<0n)neg++;else zero++;return{degreeS:A.length-1,degreeT:Math.max(...A.map(r=>r.length-1)),positive:pos,negative:neg,zero};}
function biBernsteinS(A,N=A.length-1){const choose=(n,k)=>{let v=1n;for(let i=1;i<=k;i++)v=v*BigInt(n-k+i)/BigInt(i);return v;},C=Array.from({length:N+1},()=>[0n]);for(let j=0;j<=N;j++)for(let k=0;k<=Math.min(j,A.length-1);k++)C[j]=rpAdd(C[j],A[k].map(v=>v*choose(N-k,j-k)));return C;}
function biStripPositiveT(A,candidates){let C=A.map(r=>r.slice()),factors=[];for(const f of candidates){if(f.length===1||f.some(v=>v<0n))continue;while(true){const quot=C.map(r=>rpTryDiv(r,f));if(quot.some(r=>r===null))break;C=quot;factors.push(f);}}let content=0n;for(const r of C)for(const v of r)content=riGcd(content,v);if(content>1n)C=C.map(r=>r.map(v=>v/content));return{core:biTrim(C),factors,content};}

const certified=(()=>{

const certs=[],order=[1,2,4,0,3];let positive=0,negative=0;
for(let fi=0;fi<3;fi++)for(let gi=0;gi<3;gi++){
 const data=sfRFUniversalMeanPoint(fi,gi),Q=order.map(i=>order.map(j=>data.P[i][j])),minors=[];
 for(let n=1;n<=5;n++){const poly=biDet(Q.slice(0,n).map(r=>r.slice(0,n))),stats=biStats(poly);if(stats.negative||!poly[0].some(v=>v>0n)||poly[0].some(v=>v<0n))throw Error("universal positive coefficient certificate failed");
 positive+=stats.positive;negative+=stats.negative;minors.push({n,stats,poly});}
 certs.push({...data,minors});
}
return JSON.parse(JSON.stringify({certs,order,positive,negative},(_,v)=>typeof v==="bigint"?v.toString():v));

})();
const wrong=sfRFUniversalMeanPoint(0,0,sfq(2)),wrongStats=biStats(biDet(wrong.P));if(!wrongStats.negative)throw Error('wrong center control did not fail');
return{...certified,wrongCenter:{negativeCoefficients:wrongStats.negative,expectedFailure:true}};
})();
const verified=(()=>{

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


function sfQTriangleLocalD(wt,x,mu){
 const t=wt.sq(),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),avg=v=>sfdot(mu,v),F=avg(x.map(s=>sfq(1).div(K.add(s)))),B=avg(x.map(s=>sfq(1).div(k.add(s)))),F2=avg(x.map(s=>sfq(1).div(K.add(s).sq()))),B2=avg(x.map(s=>sfq(1).div(k.add(s).sq()))),a=sfq(1).div(sfq(2).add(F)),gam=sfq(1).div(sfq(2).add(t.div(delta)).add(B)),vf=F2.sub(F.sq()),vg=B2.sub(B.sq()),cov=B.sub(F).div(K.sub(k)).sub(F.mul(B));
 const D=x.map(s=>{const d=K.add(s),e=k.add(s),L=[t.div(d),t.div(d),wt.div(e.mul(2))],aa=t.mul(s.add(sfq(1,2))).div(d),ac=t.mul(a).mul(s).mul(2).div(d),az=K.mul(gam).mul(wt).mul(s).mul(2).div(delta.mul(e)),cz=wt.mul(s.add(sfq(1,2))).div(e.mul(2)),Al=[[aa,ac,az],[ac,aa,cz],[az,cz,sfq(1)]],qa=t.mul(a).mul(s).div(d).sq().mul(vf),qz=gam.mul(wt).mul(s).div(e.mul(2)).sq().mul(vg),qcz=t.mul(a).mul(gam).mul(wt).mul(s.sq()).mul(cov).div(d.mul(e).mul(2)),O=[[qa,sfq(0),sfq(0)],[sfq(0),qa,qcz],[sfq(0),qcz,qz]];
 return Al.map((r,i)=>r.map((v,j)=>v.sub(s.add(sfq(1,4)).mul(L[i]).mul(L[j])).sub(O[i][j].mul(8))));
 });
 return{D,F,F2,B,B2,a,gam,K,k,delta,vf,vg,cov};
}
function sfQTriangleCollective12(wt,stats,M){
 const [Af,Af2,Ag,p,Cf,Cf2,Cg,z,Zf,Zg,Zg2,rho]=M,t=wt.sq(),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),diff=K.sub(k),{F,F2,B,B2}=stats,a=sfq(1).div(sfq(2).add(F)),gam=sfq(1).div(sfq(2).add(t.div(delta)).add(B)),vf=F2.sub(F.sq()),vg=B2.sub(B.sq()),fg=B.sub(F).div(diff),cov=fg.sub(F.mul(B)),bet=wt.mul(rho).add(z),b=t.mul(rho).add(wt.mul(z)),J=wt.mul(rho).div(delta.mul(2)).add(bet.mul(B).div(2)).add(z).sub(t.mul(Zg).div(2)),LX=t.mul(p).neg().div(K),LY=K.mul(2).mul(rho).add(wt.mul(z)).neg().div(t.add(1));
 const Acol=t.sq().mul(a).mul(Af.sq().add(Cf.sq())).sub(t.sq().mul(p.sq()).div(K.sq().mul(2))).add(t.mul(a).mul(p).mul(Af).mul(2)).add(K.mul(gam).mul(wt).mul(bet).mul(Ag).mul(2).div(delta)).sub(p.mul(rho).div(delta.mul(2))).add(wt.mul(bet).mul(Cg).div(2)).add(gam.mul(wt).mul(p.div(delta).sub(Cg)).mul(J)).add(rho.sq()).add(p.mul(rho)).sub(LX.mul(LY)).sub(LY.sq().div(2));
 const wf2=v1=>v1[0].sub(K.sub(sfq(1,4)).mul(v1[1])),wfg=(v1,v2)=>K.sub(sfq(1,4)).mul(v1).sub(k.sub(sfq(1,4)).mul(v2)).div(diff),bf=t.mul(p).neg().div(K.mul(2)),bg=b.neg().div(k.mul(4)),Xf=t.mul(wf2([Af,Af2]).add(wf2([Cf,Cf2]))).add(wt.mul(wfg(Zf,Zg)).div(2)),Xg=t.mul(wfg(Af,Ag).add(wfg(Cf,Cg))).add(wt.mul(Zg.sub(k.sub(sfq(1,4)).mul(Zg2))).div(2)),EV=t.mul(Af.add(Cf)).add(wt.mul(Zg).div(2)).add(bf.mul(F)).add(bg.mul(B));
 const Bcol=bf.mul(Xf).mul(2).add(bg.mul(Xg).mul(2)).add(bf.sq().mul(wf2([F,F2]))).add(bf.mul(bg).mul(wfg(F,B)).mul(2)).add(bg.sq().mul(B.sub(k.sub(sfq(1,4)).mul(B2)))).add(EV.sq().div(4));
 const ua=K.mul(Af).neg(),uc=p.sub(K.mul(Cf)),q2=uc.add(p.mul(F).div(2)),h=wt.mul(z.sub(k.mul(Zg))).div(2).add(b.mul(B).div(4)),a2=t.mul(a).sq(),covQA=Af.sub(K.mul(Af2)).sub(F.mul(ua)),covP2=Cf.sub(K.mul(Cf2)).add(p.mul(F2).div(2)).sub(F.mul(q2)),sZfg=K.mul(Zf).sub(k.mul(Zg)).div(diff),sCfg=K.mul(Cf).sub(k.mul(Cg)).div(diff),covP3=wt.mul(sZfg).div(2).add(b.mul(fg).div(4)).sub(F.mul(h)),covQAg=K.mul(Af).sub(k.mul(Ag)).div(diff).sub(B.mul(ua));
 const Rcol=a2.mul(vf).mul(ua.sq()).neg().add(a2.mul(vf).mul(p.mul(Cf.sub(K.mul(Cf2))).add(p.sq().mul(F2).div(4)).sub(q2.sq()))).add(gam.sq().mul(vg).mul(wt.mul(b).mul(Zg.sub(k.mul(Zg2))).div(4).add(b.sq().mul(B2).div(16)).sub(h.sq()))).add(a2.mul(covP2).mul(covQA).mul(2)).add(t.mul(a).mul(gam).mul(covP3).mul(covQAg).mul(2)).add(t.mul(a).mul(gam).mul(cov).mul(b.mul(sCfg).div(4).add(wt.mul(p).mul(sZfg).div(4)).add(p.mul(b).mul(fg).div(8)).sub(q2.mul(h))).mul(2));
 return{Acol,Bcol,Rcol,Ecol:Acol.sub(Bcol).sub(Rcol.mul(8))};
}
function sfQTriangle12Moments(wt,x,mu,A,C,Z,rho){
 const K=wt.sq().add(sfq(1,2)),k=wt.sq().add(1).div(2),avg=v=>sfdot(mu,v);
 const mf=(v,p)=>avg(v.map((a,i)=>a.div(p===2?K.add(x[i]).sq():K.add(x[i])))),mg=(v,p)=>avg(v.map((a,i)=>a.div(p===2?k.add(x[i]).sq():k.add(x[i]))));
 return[mf(A,1),mf(A,2),mg(A,1),avg(C),mf(C,1),mf(C,2),mg(C,1),avg(Z),mf(Z,1),mg(Z,1),mg(Z,2),rho];
}

const sfQZero=(n,m=n)=>Array.from({length:n},()=>Array(m).fill(sfq(0)));
const sfQAdd=(A,B)=>A.map((r,i)=>r.map((v,j)=>v.add(B[i][j])));
const sfQSub=(A,B)=>A.map((r,i)=>r.map((v,j)=>v.sub(B[i][j])));
const sfQScale=(A,c)=>A.map(r=>r.map(v=>v.mul(c)));
const sfQDiag=d=>d.map((v,i)=>d.map((_,j)=>i===j?v:sfq(0)));
const sfQCong=(Q,A)=>sfQM(sfQT(Q),sfQM(A,Q));
function sfQMeanSchurReduction(wt,x,mu){
 const t=wt.sq(),eps=wt,K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),avg=v=>sfdot(mu,v),n=x.length,nf=3*(n-1),dim=nf+3;
 if(!sfsum(mu).eq(1)||mu.some(p=>!p.gt(0)))throw Error("probability");
 const loc=sfQTriangleLocalD(wt,x,mu),{F,F2,B,B2,a:alpha,gam,vf,vg}=loc,stats={F,F2,B,B2},d=x.map(s=>K.add(s)),e=x.map(s=>k.add(s)),f=d.map(v=>sfq(1).div(v)),g=e.map(v=>sfq(1).div(v)),Ed=avg(d),Ed2=avg(d.map(v=>v.sq())),cf=Ed.div(Ed2),aF2=avg(d.map(v=>sfq(1).sub(cf.mul(v)).sq())),cd=avg(d.map((v,i)=>v.mul(g[i]))).div(Ed),bD2=avg(d.map((v,i)=>v.mul(g[i].sub(cd).sq()))),h=f.map((v,i)=>gam.mul(g[i]).sub(v)),cl=avg(d.map((v,i)=>v.sq().mul(h[i]))).div(Ed2),bL2=avg(d.map((v,i)=>v.sq().mul(h[i].sub(cl).sq()))),u=K.mul(gam).mul(2).div(delta).sub(sfq(1,2)),ra=sfq(1).sub(alpha.mul(2));
 const sharp=x.map((s,i)=>{
 const ds=d[i],es=e[i],w=s.add(sfq(1,4)),cp=u.mul(s).add(w.mul(s.add(sfq(1,2))).div(ds)).div(es),cm=u.mul(s).sub(sfq(1,4)).div(es),LL=s.add(sfq(1,2)).add(alpha.mul(s).mul(2)).mul(2).div(ds).sub(t.mul(w).mul(4).div(ds.sq())).sub(t.mul(alpha.sq()).mul(vf).mul(s.sq()).mul(64).div(ds.sq())).sub(t.mul(ra).mul(aF2).add(sfq(1,16)).div(ds.sq())),DD=s.add(sfq(1,2)).sub(alpha.mul(s).mul(2)).mul(2).div(ds).sub(t.mul(alpha.sq()).mul(vf).mul(s.sq()).mul(32).div(ds.sq())).sub(sfq(1).div(ds.mul(16))),ZZ=sfq(1).sub(t.mul(w).div(es.sq().mul(4))).sub(t.mul(gam.sq()).mul(vg).mul(s.sq()).mul(4).div(es.sq())).sub(t.mul(vg).div(16)).sub(t.sq().mul(vg).mul(bL2.add(gam.sq().mul(bD2))));
 return[[LL,sfq(0),cp],[sfq(0),DD,cm],[cp,cm,ZZ]];
 });
 const Psi=[[sfq(1).div(eps),sfq(1).div(eps),sfq(0)],[sfq(1).div(eps),sfq(-1).div(eps),sfq(0)],[sfq(0),sfq(0),sfq(1)]],Xi=[[sfq(0),sfq(0),sfq(0)],[sfq(1),sfq(0),sfq(0)],[sfq(0),sfq(1),sfq(0)]],Phi=x.map((_,i)=>[[f[i],sfq(0),sfq(0)],[f[i].sq(),sfq(0),sfq(0)],[g[i],sfq(0),sfq(0)],[sfq(0),sfq(1),sfq(0)],[sfq(0),f[i],sfq(0)],[sfq(0),f[i].sq(),sfq(0)],[sfq(0),g[i],sfq(0)],[sfq(0),sfq(0),sfq(1)],[sfq(0),sfq(0),f[i]],[sfq(0),sfq(0),g[i]],[sfq(0),sfq(0),g[i].sq()],[sfq(0),sfq(0),sfq(0)]]);
 const basis=Array.from({length:12},(_,i)=>Array.from({length:12},(_,j)=>sfq(i===j?1:0))),form=v=>sfQTriangleCollective12(wt,stats,v).Ecol,val=basis.map(form),J=sfQDiag(val);
 for(let i=0;i<12;i++)for(let j=i+1;j<12;j++)J[i][j]=J[j][i]=form(basis[i].map((v,k)=>v.add(basis[j][k]))).sub(val[i]).sub(val[j]).div(2);
 let R=sfQZero(12,3);for(let i=0;i<n;i++)R=sfQAdd(R,sfQScale(sfQM(Phi[i],Xi),mu[i]));R[11][2]=sfq(1);
 let C=sfQAdd(sfQCong(R,J),sfQDiag([sfq(0),sfq(0),sfq(1)]));
 const W=[],Inv=[],QQ=[],V=[],GG=sfQZero(12,dim);let total=sfQZero(dim),flow=sfQZero(nf),cross=sfQZero(nf,3),N=sfQZero(3),T=sfQZero(3),U=sfQZero(3);
 for(let i=0;i<n;i++){
 C=sfQAdd(C,sfQScale(sfQCong(Xi,loc.D[i]),mu[i]));
 W[i]=sfQAdd(sfQM(sfQT(Psi),sfQM(loc.D[i],Xi)),sfQM(sfQT(sfQM(Phi[i],Psi)),sfQM(J,R)));
 Inv[i]=sfQInv(sharp[i]);
 N=sfQAdd(N,sfQScale(Inv[i],mu[i]));T=sfQAdd(T,sfQScale(sfQM(Inv[i],W[i]),mu[i]));U=sfQAdd(U,sfQScale(sfQCong(W[i],Inv[i]),mu[i]));
 const Q=sfQZero(3,nf);for(let a=0;a<3;a++)for(let j=0;j<n-1;j++)Q[a][a*(n-1)+j]=i===j?sfq(1):i===n-1?mu[j].neg().div(mu[n-1]):sfq(0);QQ.push(Q);
 flow=sfQAdd(flow,sfQScale(sfQCong(Q,sharp[i]),mu[i]));cross=sfQAdd(cross,sfQScale(sfQM(sfQT(Q),W[i]),mu[i]));
 const Vi=sfQM(Psi,Q).map((r,a)=>[...r,...Xi[a]]);V.push(Vi);
 total=sfQAdd(total,sfQScale(sfQCong(Vi,loc.D[i]),mu[i]));const Gi=sfQM(Phi[i],Vi);for(let a=0;a<12;a++)for(let b=0;b<dim;b++)GG[a][b]=GG[a][b].add(mu[i].mul(Gi[a][b]));
 }
 GG[11][nf+2]=sfq(1);total=sfQAdd(total,sfQCong(GG,J));total[nf+2][nf+2]=total[nf+2][nf+2].add(1);
 const S=sfQAdd(sfQSub(C,U),sfQCong(T,sfQInv(N))),finiteS=sfQSub(C,sfQCong(cross,sfQInv(flow))),Factual=total.slice(0,nf).map(r=>r.slice(0,nf)),gap=sfQSub(Factual,flow);
 return {t,K,loc,stats,sharp,Inv,N,T,U,C,W,S,finiteS,flow,cross,total,Factual,gap,V,QQ,nf,dim,coefficientBounds:{aF2,bD2,bL2},original:v=>{
 const A=V.map(Vi=>sfdot(Vi[0],v)),CC=V.map(Vi=>sfdot(Vi[1],v)),Z=V.map(Vi=>sfdot(Vi[2],v));
 return sfQTriangleRemainder([wt],x,mu,A,CC,[Z],[v[nf+2]]).remainder;
 }};
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

function sfQUniversalEnergyGuard(wt,x,mu){
 const t=wt.sq(),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),F=sfdot(mu,x.map(s=>sfq(1).div(K.add(s)))),B=sfdot(mu,x.map(s=>sfq(1).div(k.add(s)))),alpha=sfq(1).div(sfq(2).add(F)),gam=sfq(1).div(sfq(2).add(t.div(delta)).add(B)),P=x.map(s=>sfQUniversalMeanPoint(wt,s,F,alpha,B,gam).P),u=x.length,nf=3*(u-1),n=nf+2;
 const field=v=>[...v,sfdot(mu.slice(0,u-1),v).neg().div(mu[u-1])],form=v=>{const L=field(v.slice(0,u-1)),D=field(v.slice(u-1,2*(u-1))),Z0=field(v.slice(2*(u-1),nf)),p=v[nf],q=v[nf+1],z=wt.mul(p).neg().div(K.mul(2)).sub(wt.mul(q).div(t.add(1))),rho=t.mul(p).neg().div(K.mul(2)).add(q.div(t.add(1))),A=L.map((a,i)=>a.add(D[i]).div(wt)),C=L.map((a,i)=>a.sub(D[i]).div(wt).add(p)),Z=Z0.map(a=>a.add(z)),E=sfQTriangleRemainder([wt],x,mu,A,C,[Z],[rho]).remainder,lower=sfdot(mu,x.map((_,i)=>{const w=[L[i],D[i],Z0[i],p,q];return sfdot(w,sfMatVec(P[i],w));}));return{E,lower};};
 const basis=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>sfq(i===j?1:0))),values=basis.map(form),mat=which=>{const M=sfQDiag(values.map(v=>v[which]));for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){const v=form(basis[i].map((v,k)=>v.add(basis[j][k])));M[i][j]=M[j][i]=v[which].sub(values[i][which]).sub(values[j][which]).div(2);}return M;},M=mat("E"),lower=mat("lower"),gap=sfQSub(M,lower),I=sfQInertia(gap),IL=sfQInertia(lower);
 if(I.negative||IL.negative||IL.zero)throw Error("universal point bound against original energy");
 return{sigma:t.str(),u,dimension:n,gap:I,lower:IL};
}
function sfQBarycentric(vertices,x,y){
 const A=vertices[0],B=vertices[1],C=vertices[2],dx=x.sub(C.x),dy=y.sub(C.y),a=A.x.sub(C.x),b=B.x.sub(C.x),c=A.y.sub(C.y),d=B.y.sub(C.y),det=a.mul(d).sub(b.mul(c)),w0=dx.mul(d).sub(b.mul(dy)).div(det),w1=a.mul(dy).sub(dx.mul(c)).div(det);return[w0,w1,sfq(1).sub(w0).sub(w1)];
}
function sfQUniversalConcavityGuard(wt,x,mu,s){
 const t=wt.sq(),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),F=sfdot(mu,x.map(v=>sfq(1).div(K.add(v)))),B=sfdot(mu,x.map(v=>sfq(1).div(k.add(v)))),alpha=sfq(1).div(sfq(2).add(F)),gam=sfq(1).div(sfq(2).add(t.div(delta)).add(B)),TA=sfQTriangleInverse(sfq(1).div(K.add(1)),sfq(1).div(K),sfq(2)),TG=sfQTriangleInverse(sfq(1).div(k.add(1)),sfq(1).div(k),sfq(2).add(t.div(delta))),wa=sfQBarycentric(TA,F,alpha),wg=sfQBarycentric(TG,B,gam),P=sfQUniversalMeanPoint(wt,s,F,alpha,B,gam).P;
 if(wa.some(v=>!v.ge(0))||wg.some(v=>!v.ge(0)))throw Error("inverse triangle enclosure");
 let mix=sfQZero(5);for(let fi=0;fi<3;fi++)for(let gi=0;gi<3;gi++)mix=sfQAdd(mix,sfQScale(sfQUniversalMeanPoint(wt,s,TA[fi].x,TA[fi].y,TG[gi].x,TG[gi].y).P,wa[fi].mul(wg[gi])));
 const gap=sfQSub(P,mix),I=sfQInertia(gap);if(I.negative)throw Error("joint matrix concavity");
 return{sigma:t.str(),s:s.str(),barycentricWeights:wa.map(v=>v.str()).concat(wg.map(v=>v.str())),gap:I};
}

const assembly=(()=>{

let entries=0,determinants=0;const profile=[];
for(const data of CERTS.certs){
 for(const wt of[sfq(1,2),sfq(1),sfq(3)])for(const s of[sfq(0),sfq(1,5),sfq(1)]){
 const t=wt.sq(),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),A=sfQTriangleInverse(sfq(1).div(K.add(1)),sfq(1).div(K),sfq(2))[data.fi],G=sfQTriangleInverse(sfq(1).div(k.add(1)),sfq(1).div(k),sfq(2).add(t.div(delta)))[data.gi],exact=sfQUniversalMeanPoint(wt,s,A.x,A.y,G.x,G.y).scaled,den=sfQEvalT(data.Den,t);
 for(let i=0;i<5;i++)for(let j=0;j<5;j++){if(!sfQEvalTS(data.P[i][j],t,s).div(den).eq(exact[i][j]))throw Error("independent polynomial assembly");entries++;}
 const Q=CERTS.order.map(i=>CERTS.order.map(j=>exact[i][j]));for(const minor of data.minors){let dn=sfq(1);for(let h=0;h<minor.n;h++)dn=dn.mul(den);const d=sfdet(Q.slice(0,minor.n).map(r=>r.slice(0,minor.n))),fromPoly=sfQEvalTS(minor.poly,t,s).div(dn);if(!d.eq(fromPoly)||!d.gt(0))throw Error("independent leading-minor identity");determinants++;}
 }
 profile.push({fi:data.fi,gi:data.gi,allFiveMinorsPositive:true});
}
return{status:"PASS",entries,determinants,profile};

})();
const energy=(()=>{

const out=[],concavity=[];
for(const wt of[sfq(1,10),sfq(1,2),sfq(1),sfq(3),sfq(10)]){
 for(const profile of[0,1]){const x=profile===0?[sfq(0),sfq(1,5),sfq(1)]:[sfq(1,6),sfq(1,3),sfq(1,2)],mu=profile===0?[sfq(1,5),sfq(3,10),sfq(1,2)]:x;out.push(sfQUniversalEnergyGuard(wt,x,mu));concavity.push(sfQUniversalConcavityGuard(wt,x,mu,sfq(2,5)));}
}return{status:"PASS",energyProfiles:out,concavityProfiles:concavity};
})();
const wrong=(()=>{

const wt=sfq(100),t=wt.sq(),s=sfq(0),K=t.add(sfq(1,2)),k=t.add(1).div(2),delta=t.mul(3).add(1).div(2),A=sfQTriangleInverse(sfq(1).div(K.add(1)),sfq(1).div(K),sfq(2))[0],G=sfQTriangleInverse(sfq(1).div(k.add(1)),sfq(1).div(k),sfq(2).add(t.div(delta)))[0],bad=sfQUniversalMeanPoint(wt,s,A.x,A.y,G.x,G.y).P,extra=wt.mul(t).div(t.add(1).mul(t.add(1)).mul(t.add(1)).mul(4));bad[0][3]=bad[3][0]=bad[0][3].add(extra);const inertia=sfQInertia(bad);if(!inertia.negative)throw Error("wrong center exact point control");return{status:"PASS",sigma:t.str(),s:s.str(),expectedIndefinite:inertia};

})();
return{assembly,energy,wrong};
})();
return{status:'PASS',theorem:'For every probability measure on [0,1] and every sigma>0, the full m1 energy A-B-8R is strictly positive on the constrained nonzero amplitude space.',scope:'Computer-assisted exact finite algebra proof for m1. Multiple letter-2 copies, full SF, the remaining middle conditions and unrestricted all-copy NPT remain open.',certificate:{families:CERTS.certs.length,leadingMinors:45,positiveCoefficients:CERTS.positive,negativeCoefficients:CERTS.negative,denominatorDegrees:CERTS.certs.map(c=>c.denDegree),wrongCenter:CERTS.wrongCenter},independent:verified};
