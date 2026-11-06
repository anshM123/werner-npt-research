function verifyRankPreservingRealification(){
function verifyPhaseBellRealificationGram(){
 const z=()=>[0n,0n],wp=[[1n,0n],[0n,1n],[-1n,-1n]];
 const add=(a,b)=>[a[0]+b[0],a[1]+b[1]],neg=a=>[-a[0],-a[1]];
 const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]-a[1]*b[1]];
 const conj=a=>[a[0]-a[1],-a[1]],sc=(a,c)=>[a[0]*c,a[1]*c],eq=(a,b)=>a[0]===b[0]&&a[1]===b[1];
 const f=Array.from({length:9},z),g=Array.from({length:9},z);
 for(let i=0;i<3;i++){f[4*i]=wp[i];g[4*i]=conj(wp[i]);}
 function inner(a,b){let o=z();for(let i=0;i<a.length;i++)o=add(o,mul(conj(a[i]),b[i]));return o;}
 if(!eq(inner(f,f),[3n,0n])||!eq(inner(g,g),[3n,0n])||!eq(inner(f,g),z()))throw Error("Bell pair not orthogonal");
 function outer(a,b){return a.flatMap(x=>b.map(y=>mul(x,conj(y))));}
 const A=[outer(f,f),outer(g,g),outer(f,g),outer(g,f)];
 function twiceM(X,stride){
  const out=X.map(x=>sc(x,2n));
  for(let i=0;i<9;i++)if(Math.floor(i/stride)%3===0)for(let j=0;j<9;j++)if(Math.floor(j/stride)%3===0){
   let tr=z();for(let k=0;k<3;k++)tr=add(tr,X[(i+k*stride)*9+j+k*stride]);
   for(let k=0;k<3;k++)out[(i+k*stride)*9+j+k*stride]=add(out[(i+k*stride)*9+j+k*stride],neg(tr));
  }
  return out;
 }
 const M=A.map(x=>twiceM(twiceM(x,3),1)),G=A.map(x=>M.map(y=>inner(x,y)));
 const E=[[33n,-3n,0n,0n],[-3n,33n,0n,0n],[0n,0n,24n,0n],[0n,0n,0n,24n]];
 for(let i=0;i<4;i++)for(let j=0;j<4;j++)if(!eq(G[i][j],[E[i][j],0n]))throw Error("Flag Gram mismatch");
 return{arithmetic:"exact BigInt in Z[omega], omega^2+omega+1=0",normalizedFlagVectors:"f=(|00>+omega|11>+omega^2|22>)/sqrt(3), g=conj(f), <f,g>=0",basis:"|f><f|, |g><g|, |f><g|, |g><f|",gramTimes36:G.map(r=>r.map(x=>x[0].toString())),constants:{a:"11/12",b:"-1/12",c:"2/3"},errorRatios:["1/11","8/11"]};
}
 const flag=verifyPhaseBellRealificationGram();
 const z=()=>[0n,0n],wp=[[1n,0n],[0n,1n],[-1n,-1n]];
 const add=(a,b)=>[a[0]+b[0],a[1]+b[1]],sc=(a,c)=>[a[0]*c,a[1]*c];
 const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]-a[1]*b[1]];
 const conj=a=>[a[0]-a[1],-a[1]],eq=(a,b)=>a[0]===b[0]&&a[1]===b[1];
 function inner(a,b){let out=z();for(let i=0;i<a.length;i++)out=add(out,mul(conj(a[i]),b[i]));return out;}
 function outer(a,b,transpose=false){return a.flatMap(x=>b.map(y=>mul(x,transpose?y:conj(y))));}
 function sumOuter(A,B,tr=false){const out=Array.from({length:A[0].length**2},z);for(let s=0;s<A.length;s++){const C=outer(A[s],B[s],tr);for(let i=0;i<out.length;i++)out[i]=add(out[i],C[i]);}return out;}
 function twoPowerM(X,n){
  const dim=3**n;let Y=X;
  for(let s=3**(n-1);s>=1;s/=3){
   const out=Y.map(x=>sc(x,2n));
   for(let i=0;i<dim;i++)if(Math.floor(i/s)%3===0)for(let j=0;j<dim;j++)if(Math.floor(j/s)%3===0){
    let tr=z();for(let k=0;k<3;k++)tr=add(tr,Y[(i+k*s)*dim+j+k*s]);
    for(let k=0;k<3;k++){const t=(i+k*s)*dim+j+k*s;out[t]=add(out[t],sc(tr,-1n));}
   }Y=out;
  }return Y;
 }
 const f=Array.from({length:9},z),g=Array.from({length:9},z);
 for(let i=0;i<3;i++){f[4*i]=wp[i];g[4*i]=conj(wp[i]);}
 function realColumn(a){
  const out=a.flatMap(x=>f.map((y,i)=>add(mul(x,y),mul(conj(x),g[i]))));
  if(out.some(x=>x[1]!==0n))throw Error("Column not real");return out;
 }
 const F=Array.from({length:3},(_,j)=>Array.from({length:3},(_,i)=>wp[(i*j)%3]));
 const I=Array.from({length:3},(_,j)=>Array.from({length:3},(_,i)=>[i===j?1n:0n,0n]));
 const cases=[{name:"complex non-Hermitian rank 2",A:[F[0],F[1]],B:[F[0],F[2]],norm:3n},{name:"complex PSD rank 2",A:[F[0],F[1]],B:[F[0],F[1]],norm:3n},{name:"negative rank-3 control only",A:I,B:I,norm:1n}];
 const records=[];
 for(const ca of cases){
  const {A,B}=ca,C=sumOuter(A,B),T=sumOuter(A,B,true),RA=A.map(realColumn),RB=B.map(realColumn),R=sumOuter(RA,RB);
  for(const cs of [A,B])for(let i=0;i<cs.length;i++)for(let j=0;j<cs.length;j++)if(!eq(inner(cs[i],cs[j]),[i===j?ca.norm:0n,0n]))throw Error("Input columns not orthogonal");
  for(const cs of [RA,RB])for(let i=0;i<cs.length;i++)for(let j=0;j<cs.length;j++)if(!eq(inner(cs[i],cs[j]),[i===j?6n*ca.norm:0n,0n]))throw Error("Rank not preserved");
  const qc=inner(C,twoPowerM(C,1)),qt=inner(T,twoPowerM(T,1)),bc=inner(C,twoPowerM(C.map(conj),1)),qr=inner(R,twoPowerM(R,3));
  if(qc[1]||qt[1]||qr[1])throw Error("Non-real q");
  const expected=66n*qc[0]-3n*(2n*bc[0]-bc[1])+48n*qt[0];
  if(qr[0]!==expected)throw Error("Independent dilation identity failed");
  if(!eq(inner(R,R),sc(inner(C,C),36n)))throw Error("Frobenius norm not preserved after normalizing R/6");
  records.push({name:ca.name,rank:A.length,twiceInputQ:qc[0].toString(),eightTimesUnnormalisedOutputQ:qr[0].toString(),identity:true});
 }
 if(records[2].eightTimesUnnormalisedOutputQ!=="-324")throw Error("Negative control failed");
 // A complex balanced 3-copy plane, realified using two Bell pairs, remains balanced on 7 copies.
 const words=[[[0,0,0],[1,1,1],[2,2,2]],[[0,1,2],[1,2,0],[2,0,1]]];
 const states=words.map(ws=>ws.flatMap((w,i)=>{
  const out=[];for(let a=0;a<3;a++)for(let b=0;b<3;b++)out.push({word:[...w,a,a,b,b],c:(i+a+b)%3===0?2n:-1n});return out;
 }));
 for(let a=0;a<2;a++)for(let b=0;b<2;b++){
  let gram=0n;for(const x of states[a])for(const y of states[b])if(x.word.every((v,i)=>v===y.word[i]))gram+=x.c*y.c;
  if(gram!==(a===b?54n:0n))throw Error("Balanced example isometry failed");
  for(let s=0;s<7;s++){
   const C=Array(9).fill(0n);for(const x of states[a])for(const y of states[b])if(x.word.every((v,i)=>i===s||v===y.word[i]))C[3*x.word[s]+y.word[s]]+=x.c*y.c;
   for(let i=0;i<3;i++)for(let j=0;j<3;j++)if(C[3*i+j]!== (a===b&&i===j?18n:0n))throw Error("Realification lost balance");
  }
 }
 return{flag,independentDirectChecks:records,balancedExample:{inputCopies:3,addedCopies:4,outputCopies:7,rank:2,allOneCopyCompressions:"I_2 delta_ij/3 exactly"},scope:"All-copy reduction with added copies; no fixed-copy restriction to real matrices and no proof of endpoint positivity."};
}
