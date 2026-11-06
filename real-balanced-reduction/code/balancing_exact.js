function verifySymmetryInterpolation(){
 const mod=x=>(x%3+3)%3,zero=()=>[0n,0n],wp=[[1n,0n],[0n,1n],[-1n,-1n]];
 const add=(a,b)=>[a[0]+b[0],a[1]+b[1]],mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]-a[1]*b[1]];
 const conj=a=>[a[0]-a[1],-a[1]],sc=(a,k)=>[a[0]*k,a[1]*k],eq=(a,b)=>a[0]===b[0]&&a[1]===b[1];
 const inner=(A,B)=>A.reduce((s,x,i)=>add(s,mul(conj(x),B[i])),zero());
 const outer=(a,b)=>a.flatMap(x=>b.map(y=>mul(x,conj(y))));
 function scaledMap(C,d,stride,p,q){
  const out=Array.from({length:d*d},zero);
  for(let i=0;i<d;i++)for(let j=0;j<d;j++){
   const a=Math.floor(i/stride)%3,b=Math.floor(j/stride)%3,ii=i+(b-a)*stride,jj=j+(a-b)*stride;
   let z=add(sc(C[i*d+j],q+p),sc(C[ii*d+jj],q-p));
   if(a===b)for(let k=0;k<3;k++)z=add(z,sc(C[(i+(k-a)*stride)*d+j+(k-b)*stride],-q));
   out[i*d+j]=z;
  }return out;
 }
 function flagGram(kind,p,q){
  const f=Array.from({length:9},zero),g=Array.from({length:9},zero);
  for(let i=0;i<3;i++)if(kind==="phase"){f[4*i]=wp[i];g[4*i]=conj(wp[i]);}else{f[4*i]=[1n,0n];g[3*i+(i+1)%3]=[1n,0n];}
  if(!eq(inner(f,f),[3n,0n])||!eq(inner(g,g),[3n,0n])||!eq(inner(f,g),zero()))throw Error("Flag orthogonality failed");
  const A=[outer(f,f),outer(g,g),outer(f,g),outer(g,f)],M=A.map(X=>scaledMap(scaledMap(X,9,1,p,q),9,3,p,q));
  return A.map(X=>M.map(Y=>inner(X,Y)));
 }
 let gramChecks=0;
 for(const[p0,q0]of [[0,1],[1,1],[-1,1],[1,2],[1,3]]){
  const p=BigInt(p0),q=BigInt(q0);
  const a=3n*(4n*q*q+6n*p*q+p*p),b=3n*(4n*q*q-6n*p*q+p*p),c=12n*(q*q+p*p);
  const expectedPhase=[[a,b,0n,0n],[b,a,0n,0n],[0n,0n,c,0n],[0n,0n,0n,c]];
  const ar=3n*(7n*q*q+4n*p*p),br=-3n*q*q,cr=3n*(3n*q*q+4n*p*q+p*p),dr=3n*(3n*q*q-4n*p*q+p*p);
  const expectedReal=[[ar,br,0n,0n],[br,ar,0n,0n],[0n,0n,cr,dr],[0n,0n,dr,cr]];
  for(const[kind,E]of [["phase",expectedPhase],["real",expectedReal]]){
   const G=flagGram(kind,p,q);
   for(let i=0;i<4;i++)for(let j=0;j<4;j++){if(!eq(G[i][j],[E[i][j],0n]))throw Error("Flag Gram identity failed");gramChecks++;}
  }
 }
 // Check the one-copy complex projector and its symmetric companion in Z[omega].
 const u=Array.from({length:3},()=>[1n,0n]),v=[wp[0],wp[2],wp[1]];
 if(!eq(inner(u,u),[3n,0n])||!eq(inner(v,v),[3n,0n])||!eq(inner(u,v),zero()))throw Error("Input orthogonality failed");
 const C=outer(u,u).map((x,i)=>add(x,outer(v,v)[i]));
 const T=u.flatMap(x=>u.map(y=>mul(x,y))).map((x,i)=>add(x,v.flatMap(a=>v.map(b=>mul(a,b)))[i]));
 const barC=C.map(conj);
 for(const[p0,q0]of [[0,1],[1,2],[1,1]]){
  const p=BigInt(p0),q=BigInt(q0),MC=scaledMap(C,3,1,p,q),MT=scaledMap(T,3,1,p,q),Mb=scaledMap(barC,3,1,p,q);
  // C and T are three times their normalized operators.
  if(!eq(inner(C,MC),[9n*(p-q),0n])||!eq(inner(C,Mb),[-9n*(p+q),0n])||!eq(inner(T,MT),[27n*q,0n]))throw Error("Input Gram values failed");
 }
 // Independent residue-automaton contraction of the sparse real projector.
 const tuples=Array.from({length:81},(_,x)=>{const a=[];for(let k=0;k<4;k++){a.push(x%3);x=Math.floor(x/3);}return a;});
 const enc=a=>a.reduce((s,x,k)=>s+x*3**k,0),cos=x=>mod(x)===0?2n:-1n;
 function actualNumerator(R,p,q){
  const kernel=([i,j,k,l])=>(i===k&&j===l?q+p:0n)+(i===l&&j===k?q-p:0n)-(i===j&&k===l?q:0n);
  const weights=tuples.map(x=>kernel(x)**2n);let mass=Array(81).fill(0n);mass[0]=1n;
  for(let r=0;r<R;r++){const next=Array(81).fill(0n);for(let a=0;a<81;a++)if(mass[a])for(let b=0;b<81;b++)if(weights[b]){const t=enc(tuples[a].map((x,k)=>mod(x+tuples[b][k])));next[t]+=mass[a]*weights[b];}mass=next;}
  let total=0n;
  for(let a=0;a<81;a++)if(mass[a])for(const t of tuples){
   const K=kernel(t);if(!K)continue;const s=tuples[a];
   const Cij=cos(s[0])*cos(s[1])+cos(s[0]-t[0])*cos(s[1]-t[1]);
   const Ckl=cos(s[2])*cos(s[3])+cos(s[2]-t[2])*cos(s[3]-t[3]);
   total+=mass[a]*K*Cij*Ckl;
  }return total;
 }
 let projectorNormChecks=0;
 for(let R=1;R<=12;R++){let xx=0n,yy=0n,xy=0n;const multiplicity=3n**BigInt(R-1);for(let t=0;t<3;t++)for(let i=0;i<3;i++){xx+=multiplicity*cos(t)**2n;yy+=multiplicity*cos(t-i)**2n;xy+=multiplicity*cos(t)*cos(t-i);}const D=2n*3n**BigInt(R+1);if(xx!==D||yy!==D||xy!==0n)throw Error("Real projector factors not orthogonal");projectorNormChecks++;}
 let transferChecks=0;
 for(let R=1;R<=8;R++)for(const[p0,q0]of [[0,1],[1,3],[1,2],[1,1]]){
  const p=BigInt(p0),q=BigInt(q0),N=2*R+1,D=2n*3n**BigInt(R+1);
  const actual=actualNumerator(R,p,q),actualDen=(2n*q)**BigInt(N)*D*D;
  const A=4n*q*q+6n*p*q+p*p,B=4n*q*q-6n*p*q+p*p,Cc=4n*(q*q+p*p);
  const form=(p-q)*A**BigInt(R)-(p+q)*B**BigInt(R)+3n*q*Cc**BigInt(R),formDen=4n*q*(12n*q*q)**BigInt(R);
  if(actual*formDen!==form*actualDen)throw Error("Independent residue transfer disagrees");transferChecks++;
 }
 // Verify the real twelve-element orthogonal 1-design on every matrix unit.
 const signs=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],O=[];
 for(let t=0;t<3;t++)for(const s of signs)O.push(Array.from({length:3},(_,i)=>Array.from({length:3},(_,j)=>i===(j+t)%3?s[j]:0)));
 let designChecks=0;
 for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(let a=0;a<3;a++)for(let b=0;b<3;b++){
  const z=O.reduce((s,M)=>s+M[i][a]*M[j][b],0);if(z!==((i===j&&a===b)?4:0))throw Error("Real orthogonal 1-design failed");designChecks++;
 }
 let flagErasureChecks=0;
 for(let h=0;h<12;h++)for(let k=0;k<12;k++)if(h!==k)for(let bit=0;bit<4;bit++)for(let rep=0;rep<2;rep++)for(let physical=0;physical<2;physical++){
  let otherOrthogonalPair=false;for(let b=0;b<4;b++)for(let r=0;r<2;r++)if(!(b===bit&&r===rep)&&((h>>b)&1)!==((k>>b)&1))otherOrthogonalPair=true;
  if(!otherOrthogonalPair)throw Error("An erased flag copy preserves a cross-label branch");flagErasureChecks++;
 }
 const gcd=(a,b)=>b?gcd(b,a%b):(a<0n?-a:a),rat=(a,b)=>{const g=gcd(a,b);return(a/g).toString()+"/"+(b/g).toString();};
 const num=-(29n**5n)-3n*5n**5n+6n*20n**5n,den=8n*48n**5n;if(!(num<0n))throw Error("Concrete negative example failed");
 return{arithmetic:"exact BigInt in Z[omega] and integers",flagGramEntryChecks:gramChecks,independentResidueTransferChecks:transferChecks,projectorFactorNormChecks:projectorNormChecks,realOneDesignEntryChecks:designChecks,flagErasureCombinatorialChecks:flagErasureChecks,realProjectorFamily:"N=1+2R, P_R=(x_R x_R^T+y_R y_R^T)/(2*3^(R+1)); amplitudes 2 cos(2pi sum(k)/3) and 2 cos(2pi(sum(k)-i)/3)",closedForm:"q_(M_z)^(1+2R)(P_R)=((z-1)*a_z^R-(1+z)*b_z^R+3*c_z^R)/4",phaseConstants:{a:"(4+6z+z^2)/12",b:"(4-6z+z^2)/12",c:"(1+z^2)/3"},example:{z:"1/2",R:5,copies:11,q:rat(num,den),endpointQ:rat(3n*8n**5n+2n,4n*12n**5n)},realBalancingFlagConstants:{a:"(7+4z^2)/12",b:"-1/12",c:"(3+4z+z^2)/12",d:"(3-4z+z^2)/12",eta:"(z^2+4z+3)/(4z^2+7)<=8/11<1 for 0<=z<=1; difference numerator 21z^2-44z+23=(1-z)(23-21z)"},scope:"For every 0<z<1 some real PSD rank-two projector, and after the proved balancing construction some real balanced rank-two projector at an odd copy count, has negative M_z form. M_0 has the separate all-copy real PSD rank-two SOS; M_1 is the still-open Werner endpoint."};
}
