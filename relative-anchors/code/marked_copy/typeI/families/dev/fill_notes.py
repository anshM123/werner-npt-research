"""Development helper (not part of any proof): fills the tables of NOTES.md, Section 8, from the logs, via tables.py
and d1table.sh.  python dev/fill_notes.py   (run in the stage folder)"""
import os, re, subprocess, sys

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def pretty(row):
    row = row.replace('>=', ' ≥ ').replace('tp', 'τ′').replace('=', ' = ')
    row = re.sub(r'(?<=[\s\[(,])-(\d)', r'−\1', row)
    row = re.sub(r'(?<![\d.])(\d{1,3})(\d{3})(\d{3})(?![\d.])', r'\1 \2 \3', row)
    row = re.sub(r'(?<![\d.])(\d{1,3})(\d{3})(?![\d.])', r'\1 \2', row)
    row = re.sub(r'  +', ' ', row)
    return row


def tag(spec):
    return spec.replace('-', 'm').replace('/', '_').replace(',', '_')


def run(args):
    return subprocess.run(args, capture_output=True, text=True, encoding='utf-8').stdout.strip().split('\n')


s = open('NOTES.md', encoding='utf-8').read()
# MU Proposition D with the split
fams = open('lists/MU_sym_2d_fams.txt').read().split()
logs = [f'logs/sd_{tag(f)}.log' for f in fams if os.path.exists(f'logs/sd_{tag(f)}.log')]
rows = run([sys.executable, 'tables.py', 'sd'] + logs)
if 'TABLE_SD' in s:
    s = s.replace('TABLE_SD', '\n'.join(pretty(r) for r in rows if r))
# D2
fams = open('lists/MT_sym_3d.txt').read().split()
dirs = [f'd2/{tag(f)}' for f in fams]
rows = run([sys.executable, 'tables.py', 'd2'] + dirs)
out = []
for r in rows:
    if not r:
        continue
    m = re.match(r'\| MT/sym/([-\dn]+),([-\dn]+),n,n \|', r)
    q0, t = m.group(1), m.group(2)
    lab = '[' + ', '.join([f'q0 = {q0}' if q0 != 'n' else 'q0 ≥ 4', f't = {t}' if t != 'n' else 't ≥ 4',
                           'i ≥ 2', 'τ′ ≥ 2']) + ']'
    r = '| ' + lab + ' |' + r[m.end():]
    out.append(pretty(r).replace('q0 = −', 'q0 = −'))
if 'TABLE_D2\n' in s:
    s = s.replace('TABLE_D2\n', '\n'.join(out) + '\n', 1)
# the cheap MU families of Proposition D1 with the split
fams = open('lists/MU_sym_2e_cheap.txt').read().split()
logs = [f'logs/sd_{tag(f)}.log' for f in fams if os.path.exists(f'logs/sd_{tag(f)}.log')]
txt = ''.join(open(f, encoding='utf-8').read() for f in logs)
open('logs/sd_MU_2e_cheap_all.log', 'w', encoding='utf-8').write(txt)
rows = run(['bash', 'd1table.sh', 'logs/sd_MU_2e_cheap_all.log'])
if 'TABLE_D1MU' in s:
    head = ('| family | vectors | sliced size: 0..B | Ñ terms | degrees | negative | time (s) | result |\n'
            '|---|---|---|---|---|---|---|---|\n')
    s = s.replace('TABLE_D1MU', head + '\n'.join(pretty(r) for r in rows if r))
open('NOTES.md', 'w', encoding='utf-8').write(s)
print('filled')
