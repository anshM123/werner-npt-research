(()=>{

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


const auditCycleRatio=()=>{
 const tail=Array.from({length:12},(_,i)=>({n:i+1,ratio:null,word:null,col:null}));
 let rowRatio=null,ratioWord=null,ratioRow=null,globalMin=null,globalMax=null;
 const I=Array.from({length:9},(_,i)=>Array.from({length:9},(_,j)=>[i===j?1n:0n,0n,0n]));
 const visit=(P,w)=>{
  if(w.length===12)return;
  for(let b=0;b<3;b++){
   const Q=cycleAppend(P,b),word=w.concat(b),n=word.length;
   if(n>=8&&n<=11){
    for(let j=0;j<9;j++){let pos=0n,neg=0n;
     for(let i=0;i<9;i++)for(const v of Q[i][j]){if(v>=0n)pos+=v;else neg-=v;}
     const r=tail[n-1];if(neg>0n&&(!r.ratio||pos*BigInt(r.ratio[1])<BigInt(r.ratio[0])*neg)){r.ratio=[String(pos),String(neg)];r.word=word;r.col=j;}
    }
   }
   if(n===12){
    for(let i=0;i<9;i++){let lo=Q[i][0][0],hi=lo;
     for(let j=0;j<9;j++)for(const v of Q[i][j]){if(v<lo)lo=v;if(v>hi)hi=v;}
     if(lo<=0n)throw Error("Expected positive word");if(globalMin===null||lo<globalMin)globalMin=lo;if(globalMax===null||hi>globalMax)globalMax=hi;
     if(!rowRatio||hi*BigInt(rowRatio[1])>BigInt(rowRatio[0])*lo){rowRatio=[String(hi),String(lo)];ratioWord=word;ratioRow=i;}
    }
   }visit(Q,word);
  }
 };visit(I,[]);
 return{globalMin:String(globalMin),globalMax:String(globalMax),rowRatio12:rowRatio,ratioWord,ratioRow,tail:tail.filter(r=>r.n>=8)};
};


const cycleTailMarginAudit=()=>{
 const rows=Array.from({length:11},(_,i)=>({n:i+1,gamma:null,word:null,col:null}));
 const I=Array.from({length:9},(_,i)=>Array.from({length:9},(_,j)=>[i===j?1n:0n,0n,0n]));
 const visit=(P,w)=>{
  if(w.length===11)return;
  for(let b=0;b<3;b++){
   const Q=cycleAppend(P,b),word=w.concat(b),r=rows[word.length-1];
   if(word.length>=6)for(let j=0;j<9;j++){
    let pos=0n,neg=0n;for(let i=0;i<9;i++)for(const v of Q[i][j]){if(v>=0n)pos+=v;else neg-=v;}
    const gamma=365472n*pos-4764744n*neg;
    if(r.gamma===null||gamma<BigInt(r.gamma)){r.gamma=String(gamma);r.word=word;r.col=j;}
   }visit(Q,word);
  }
 };visit(I,[]);return rows.filter(r=>r.n>=6);
};


const r=auditCycleRatio(),g=cycleTailMarginAudit();
if(r.globalMin!=="365472")throw Error("12-block minimum");
if(BigInt(r.rowRatio12[0])*365472n!==4764744n*BigInt(r.rowRatio12[1]))throw Error("row ratio");
if(g.some(x=>BigInt(x.gamma)<=0n))throw Error("tail bound");
const G=[12,13,18,19,20,21,22,23],seen=Array(42).fill(false),D=Array(42).fill(null);seen[0]=true;D[0]=[];
for(let n=1;n<=41;n++)for(const z of G)if(n>=z&&seen[n-z]){seen[n]=true;D[n]=D[n-z].concat(z);break;}
if(!seen.slice(30).every(Boolean))throw Error("length coverage");
return{passed:true,rowAudit:r,gamma:g,coverage:D.slice(30)};

})()
