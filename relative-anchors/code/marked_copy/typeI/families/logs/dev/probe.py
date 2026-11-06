"""A measurement, not a proof: the certificate of tprove.py on one family with the pivots printed as they are taken
and a time limit; it reports how far the elimination got.  Ends with PASS when the run itself completed (whatever the
certificate status), so that it can be a step of run_extra.sh.
python probe.py FAM SECONDS"""
import sys, time
import tprove as P

case, mode, fam = P.parse_fam(sys.argv[1])
tl = float(sys.argv[2])
t0 = time.time()
res, dt, info = P.prove(case, mode, fam, tlimit=tl, verbose=True)
piv = info.get('pivots', [])
print(f'{P.fam_label(case, mode, fam)}: status {res}; {len(piv)} pivots taken of {len(info.get("names", [])) - len(info.get("gauge", []))}'
      f' in {dt:.1f}s; form {info.get("form_time")}s, largest entry {info.get("form_terms")} terms')
print('PASS')
