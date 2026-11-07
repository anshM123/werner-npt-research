"""The table of Section 1 of the notes at k = 8 from the certificate logs: per family the block, the number of states,
the known kernel, dim W, the number of configurations, the parts (and the largest part), the leading minors, the largest
numerator degree and the total time.  python k8table.py"""
import re, glob, os
from collections import OrderedDict

FAMS = ('j3m0 j3m1 j3m2 j3m3 j4m0 j4m1 j4m2 j4m3 j4m4 j5m0 j5m1 j5m2 j5m3 j5m4 j5m5 j6m0 j6m1 j6m2 j6m3 j6m4 j6m5 '
        'j6m6 j7m2 j7m3 j7m4 j7m5').split()
rows = OrderedDict()
tot_cert = tot_parts = tot_minors = 0
allpass = True
for f in FAMS:
    logs = sorted(glob.glob(f'logs/run_cert_k8_{f}_*_*.log'))
    logs = [l for l in logs if 'first_attempt' not in l]
    st = ker = dimw = None
    nconf = nparts = nminors = maxpart = maxdeg = 0
    tsum = 0.0
    for lg in logs:
        s = open(lg, encoding='utf-8', errors='replace').read()
        lines = [x.strip() for x in s.splitlines() if x.strip()]
        allpass &= bool(lines) and lines[-1] == 'PASS'
        m = re.search(r'block over Q\(q\): (\d+) states', s)
        st = int(m.group(1))
        m = re.search(r'(\d+) known kernel vectors', s)
        ker = int(m.group(1))
        m = re.search(r'W: dimension (\d+)', s)
        dimw = int(m.group(1))
        parts = [int(x) for x in re.findall(r'\(dim \d+\): size (\d+),', s)]
        nconf += 1
        nparts += len(parts)
        nminors += sum(parts)
        maxpart = max([maxpart] + parts)
        m = re.search(r'largest degrees: numerator (\d+)', s)
        maxdeg = max(maxdeg, int(m.group(1)))
        m = re.search(r'^time ([0-9.]+)s', s, re.M)
        tsum += float(m.group(1))
    j, mm = int(f[1]), int(f[3])
    block = f'(R∞)_{j}, letter-1' if mm == 0 else f'(R∞)_{j}, |T| = {mm}'
    rows[f] = (block, st, ker, dimw, nconf, nparts, maxpart, nminors, maxdeg, tsum)
    tot_cert += nconf
    tot_parts += nparts
    tot_minors += nminors
print('| family | block | states | known kernel | dim W | configurations | parts (largest) | minors | degree | time |')
print('|---|---|---|---|---|---|---|---|---|---|')
for f, (block, st, ker, dimw, nconf, nparts, maxpart, nminors, maxdeg, tsum) in rows.items():
    print(f'| {f} | {block} | {st} | {ker} | {dimw} | {nconf} | {nparts} ({maxpart}) | {nminors} | {maxdeg} | '
          f'{tsum:.0f} s |')
print()
print(f'{tot_cert} certificates, {tot_parts} parts, {tot_minors} minors; every log ends with PASS: {allpass}')
