"""A stress test of the analytic reduction of the energy theorem (NOTES.md, Section 3), using the guards of the
certificate itself: at 40 further random weight measures (20 physical, 20 formal, some with one dominant free weight;
sigma from 1/400 to 324) the full matrix inequality E >= int y^T P y dmu holds (the gap has no negative eigenvalue and
the lower form is positive definite), and the concavity/triangle step holds.  The certificate program is copied with
these profiles added to its own profile loop and run with run_js.py.   python stress_profiles.py"""
import os, subprocess, sys
src = open('one_letter_certificate.js', encoding='utf-8').read()
anchor = '}return{status:"PASS",energyProfiles:out,concavityProfiles:concavity};'
assert src.count(anchor) == 1
extra = r'''}
{let st=987654321n;const rnd=m=>{st=(st*6364136223846793005n+1442695040888963407n)%(1n<<63n);return Number(st%BigInt(m));};
const extra=[];
for(let it=0;it<40;it++){
 const u=2+rnd(5),formal=(it%2===1);
 let xs=Array.from({length:u},()=>sfq(1+rnd(400),1+rnd(20)));
 if(it%5===4)xs=xs.map((v,i)=>i===0?v.mul(sfq(1000)):v);
 let x,mu;
 if(!formal){const S=sfsum(xs);x=xs.map(v=>v.div(S));mu=x;}
 else{x=Array.from({length:u},()=>sfq(1+rnd(97),100));const ms=Array.from({length:u},()=>sfq(1+rnd(50)));const M=sfsum(ms);mu=ms.map(v=>v.div(M));}
 const wtc=[sfq(1,20),sfq(1,3),sfq(1),sfq(3),sfq(20)][rnd(5)].mul(sfq(1+rnd(9),5));
 const e=sfQUniversalEnergyGuard(wtc,x,mu);sfQUniversalConcavityGuard(wtc,x,mu,sfq(1+rnd(9),10));
 extra.push({u,formal,sigma:wtc.sq().str(),gapNegative:e.gap.negative,lowerPositiveDefinite:!e.lower.negative&&!e.lower.zero});}
out.push({extraRandomProfiles:extra.length,list:extra});}
return{status:"PASS",energyProfiles:out,concavityProfiles:concavity};'''
open('stress_certificate.js', 'w', encoding='utf-8', newline='\n').write(src.replace(anchor, extra))
r = subprocess.run([sys.executable, 'run_js.py', 'stress_certificate.js'], capture_output=True, text=True)
os.remove('stress_certificate.js')
out = r.stdout
n = out.count('"lowerPositiveDefinite": true'); g = out.count('"gapNegative": 0')
print(f'40 random measures: gap without negative eigenvalue {g}/40, lower form positive definite {n}/40; certificate status '
      f'{"PASS" if out.count(chr(34) + "status" + chr(34) + ": " + chr(34) + "PASS" + chr(34)) >= 4 else "FAIL"}')
print('PASS' if r.returncode == 0 and n == 40 and g == 40 else 'FAIL')
