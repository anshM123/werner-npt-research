(()=>{

const cycleTransferPowerAudit=(maxN)=>{
 const audits=[],certs=[];
 for(let b=0;b<3;b++){
  const edges=Array.from({length:9},()=>[]);for(let x=0;x<3;x++)for(let y=0;y<3;y++)for(let z=0;z<3;z++)edges[3*x+y].push({to:3*y+z,label:y,c:y===0&&(x+z+b)%3===0?-1n:2n});
  let P=Array.from({length:9},(_,i)=>Array.from({length:9},(_,j)=>[i===j?1n:0n,0n,0n])),rows=[];
  for(let n=1;n<=maxN;n++){
   const B=Array.from({length:9},()=>Array.from({length:9},()=>[0n,0n,0n]));for(let i=0;i<9;i++)for(let j=0;j<9;j++)for(const e of edges[j])for(let a=0;a<3;a++)B[i][e.to][(a+e.label)%3]+=P[i][j][a]*e.c;P=B;
   const vals=P.flat(2),negative=vals.filter(v=>v<0n).length,zero=vals.filter(v=>v===0n).length,min=vals.reduce((a,v)=>v<a?v:a);rows.push({n,negative,zero,min:String(min)});
  }
  let base=null;for(let m=1;2*m-1<=maxN;m++)if(rows.slice(m-1,2*m-1).every(r=>r.negative===0)){base=m;break;}
  audits.push({b,rows,base});
 }
 return audits;
};


const cycleTransferSpectra=(maxN,keep=[])=>{
 const D=9,polyz=()=>[0n,0n,0n],all=[];
 for(let b=0;b<3;b++){
  const edges=Array.from({length:D},()=>[]);
  for(let x=0;x<3;x++)for(let y=0;y<3;y++)for(let z=0;z<3;z++)edges[3*x+y].push({to:3*y+z,label:y,c:y===0&&(x+z+b)%3===0?-1n:2n});
  let P=Array.from({length:D},(_,i)=>Array.from({length:D},(_,j)=>[i===j?1n:0n,0n,0n])),vals=[];
  for(let n=1;n<=maxN;n++){
   const B=Array.from({length:D},()=>Array.from({length:D},polyz));
   for(let i=0;i<D;i++)for(let j=0;j<D;j++)for(const e of edges[j])for(let a=0;a<3;a++)B[i][e.to][(a+e.label)%3]+=P[i][j][a]*e.c;
   P=B;const z=polyz();for(let i=0;i<D;i++)for(let a=0;a<3;a++)z[a]+=P[i][i][a];vals.push(z);
  }all.push(vals);
 }
 const out=[],bad=[],saved=[];let smallestRatio=null;
 for(let n=3;n<=maxN;n++){
  const sums=Array.from({length:3},(_,a)=>Array.from({length:3},(_,b)=>all[b][n-1][a])),pairs=[[1,0],[0,1],[1,1],[1,2]],pvals=pairs.map(([a,b])=>sums[a][b]);let mn=pvals.reduce((v,x)=>v<x?v:x),z=sums[0][0],proj=4n*z+2n*mn,r1=3n*mn,r2=2n*z+mn,rank2=r1<r2?r1:r2,den=3n**BigInt(n-1)*2n**BigInt(n);
  const ratio=Number(rank2)/Number(3n*den);
  if(n>=4&&(proj<0n||rank2<0n))bad.push({n,proj:String(proj),rank2:String(rank2)});
  if(keep.includes(n))saved.push({n,sums:sums.map(r=>r.map(String)),denominator:String(den),projectorNumerator:String(proj),projectorDenominator:String(3n*den),rank2Numerator:String(rank2),rank2Denominator:String(3n*den),approxRank2:ratio});
  if(n>=4&&Number.isFinite(ratio)&&(!smallestRatio||ratio<smallestRatio.value))smallestRatio={n,value:ratio};
 }
 return{maxN,bad,smallestRatio,saved};
};


const cycleWordProduct=(word)=>{
 let P=Array.from({length:9},(_,i)=>Array.from({length:9},(_,j)=>[i===j?1n:0n,0n,0n]));
 for(const b of word)P=cycleAppend(P,b);return P;
};
const cycleAppend=(P,b)=>{
 const Q=Array.from({length:9},()=>Array.from({length:9},()=>[0n,0n,0n]));
 for(let y=0;y<3;y++)for(let z=0;z<3;z++)for(let i=0;i<9;i++)for(let a=0;a<3;a++){
  const v=P[i][y][a]+P[i][3+y][a]+P[i][6+y][a];
  Q[i][3*y+z][(a+y)%3]=2n*v-(y===0?3n*P[i][3*((6-z-b)%3)][a]:0n);
 }return Q;
};
const cycleAllWordAudit=(N)=>{
 const rows=Array.from({length:N},(_,i)=>({n:i+1,words:0,negProductWords:0,zeroProductWords:0,minCoeff:null,minTrace:null,worstWord:null}));
 const I=Array.from({length:9},(_,i)=>Array.from({length:9},(_,j)=>[i===j?1n:0n,0n,0n]));
 const visit=(P,w)=>{
  if(w.length===N)return;
  for(let b=0;b<3;b++){
   const Q=cycleAppend(P,b),word=w.concat(b),r=rows[word.length-1];r.words++;
   let m=Q[0][0][0];for(const row of Q)for(const v of row)for(const x of v)if(x<m)m=x;
   if(m<0n)r.negProductWords++;if(m===0n)r.zeroProductWords++;
   if(r.minCoeff===null||m<BigInt(r.minCoeff)){r.minCoeff=String(m);r.worstWord=word;}
   let mt=null;for(let a=0;a<3;a++){let s=0n;for(let j=0;j<9;j++)s+=Q[j][j][a];if(mt===null||s<mt)mt=s;}
   if(r.minTrace===null||mt<BigInt(r.minTrace))r.minTrace=String(mt);
   visit(Q,word);
  }
 };visit(I,[]);return rows;
};


const cycleWordSpectrum=(word)=>{
 let P=Array.from({length:9},(_,i)=>Array.from({length:9},(_,j)=>[i===j?1n:0n,0n,0n]));
 for(const b of word){
  const Q=Array.from({length:9},()=>Array.from({length:9},()=>[0n,0n,0n]));
  for(let j=0;j<9;j++){const x=(j/3)|0,y=j%3;
   for(let z=0;z<3;z++){const to=3*y+z,c=y===0&&(x+z+b)%3===0?-1n:2n;
    for(let i=0;i<9;i++)for(let a=0;a<3;a++)Q[i][to][(a+y)%3]+=c*P[i][j][a];
   }
  }P=Q;
 }return Array.from({length:3},(_,a)=>P.reduce((s,r,i)=>s+r[i][a],0n));
};


const cycleDirectCosetNumerator=(n,word,a)=>{
 const D=3**n,vec=Array.from({length:D},(_,r)=>Array.from({length:n},(_,i)=>(r/3**(n-1-i)|0)%3));
 const mod=x=>(x%3+3)%3,xx=Array(n).fill(0);xx[0]=a;
 const zz=word.map((v,i)=>mod(v+xx[(i+n-1)%n]+xx[(i+1)%n]));
 const sum=v=>v.reduce((s,x)=>s+x,0)%3,phase=v=>v.reduce((s,x,i)=>s+x*v[(i+1)%n],0)%3;
 const entries=[];
 for(let r=0;r<D;r++)for(let c=0;c<D;c++){
  const rr=vec[r],cc=vec[c],cx=cc.map((v,i)=>mod(v+xx[i]));
  if(sum(rr)!==sum(cx))continue;
  const e=mod(phase(rr)-phase(cx)+zz.reduce((s,z,i)=>s+z*cc[i],0));
  entries.push([rr,cc,e]);
 }
 let total=0n;
 const roots=[[1,0],[0,1],[-1,-1]];
 for(let mask=0;mask<1<<n;mask++){
  const J=[],K=[];for(let i=0;i<n;i++)((mask>>i)&1?J:K).push(i);
  const A=new Map();
  for(const [r,c,e]of entries){
   if(J.some(i=>r[i]!==c[i]))continue;
   let ri=0,ci=0;for(const i of K){ri=3*ri+r[i];ci=3*ci+c[i];}
   const key=ri*3**K.length+ci,v=A.get(key)||[0,0],root=roots[e];
   v[0]+=root[0];v[1]+=root[1];A.set(key,v);
  }
  let z=0n;for(const [x,y]of A.values())z+=BigInt(x*x-x*y+y*y);
  total+=BigInt((-1)**J.length*2**(n-J.length))*z;
 }return total;
};


const assert=(x,msg)=>{if(!x)throw Error(msg);};
const powers=cycleTransferPowerAudit(15);
assert(powers.every(b=>b.rows.filter(r=>r.n>=8).every(r=>r.negative===0&&r.zero===0)),"homogeneous powers");
const small=cycleTransferSpectra(7,[4,5,6,7]);
assert(small.bad.length===0,"small homogeneous traces");
const all=cycleAllWordAudit(13);
for(const n of [12,13]){
 const r=all[n-1];
 assert(r.words===3**n&&r.negProductWords===0&&r.zeroProductWords===0,"mixed-word certificate");
 assert(r.minCoeff===(n===12?"365472":"2773116"),"minimum coefficient");
}
assert(all.filter(r=>r.n>=4).every(r=>BigInt(r.minTrace)>0n),"small all-coset traces");
let seed=123;const rnd=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/2**32;};
let physical=0;
for(const n of [4,5])for(let j=0;j<(n===4?25:6);j++){
 const word=Array.from({length:n},()=>Math.floor(rnd()*3)),a=j%3;
 const z=cycleWordSpectrum(word),P=cycleWordProduct(word);
 assert(z.every((v,k)=>v===P.reduce((s,row,i)=>s+row[i][k],0n)),"independent transfer update");
 assert(cycleDirectCosetNumerator(n,word,a)===3n**BigInt(n)*z[a],"direct physical partial traces");
 physical++;
}
return {passed:true,homogeneousPowers:"8..15",allMixedWords:all.filter(r=>r.n>=12),directPhysicalChecks:physical,finiteAllCosetMinimum:all.filter(r=>r.n>=4).map(r=>({n:r.n,numerator:r.minTrace,denominator:String(3n**BigInt(r.n-1)*2n**BigInt(r.n))}))};

})()
