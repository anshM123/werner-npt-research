"""Development helper (not part of any proof): writes Section 9 of NOTES.md from logs/run_coverage_st4.log."""
import os, re
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
log = open('logs/run_coverage_st4.log', encoding='utf-8').read()


def num(x):
    x = int(x)
    s = f'{x:,}'.replace(',', ' ')
    return s


per = {}
for k, rest in re.findall(r'^k=(\d+): (.*)$', log, re.M):
    d = dict((a.strip(), int(b)) for a, b in re.findall(r'([A-Za-z0-9=\-() ]+?) (\d+)(?:,|$)', rest))
    per[int(k)] = d
tot = re.search(r'^k = 4\.\.120: (.*)$', log, re.M).group(1)
totd = dict((a.strip(), int(b)) for a, b in re.findall(r'([A-Za-z0-9=\-() ]+?) (\d+)(?:,|$)', tot))
rows = {}
for c, mo, fam, unc, sizes, kmin in re.findall(
        r'^(M[TU]) (sym|anti): 30 <= k <= 120: (\d+) forms in proved families; not covered: (\d+) \(by symbolic sizes '
        r'of the canonical family: (\{[^}]*\})\); smallest k not covered: (\S+)$', log, re.M):
    rows[(c, mo)] = (int(fam), int(unc), eval(sizes), kmin)
q3 = dict((c, (int(a), int(b))) for c, a, b in re.findall(
    r'^(M[TU]) sym: of the (\d+) forms not covered, (\d+) lie in a family', log, re.M))
fr = re.search(r'families read: (\d+) proved, (\d+) empty, (\d+) with the three-row part certified', log)
ok = log.strip().endswith('PASS')
out = []
out.append('`coverage_st3.py 120 29 30 res/*.pkl res3/q3.pkl` (`logs/run_coverage_st4.log`, which ends with '
           + ('PASS' if ok else 'FAIL') + f'), with every result of stages 3 and 4 ({fr.group(1)} families certified as a '
           f'whole; the {fr.group(3)} three-row certificates of Section 8.6 cover no form and are only counted):')
out.append('')
out.append('| | 30 ≤ k ≤ 120, in certified families | 30 ≤ k ≤ 120, open | open, by symbolic sizes of the canonical family | '
           'open forms whose Q₃ is certified |')
out.append('|---|---|---|---|---|')
for c in ('MT', 'MU'):
    for mo in ('sym', 'anti'):
        fam, unc, sizes, kmin = rows[(c, mo)]
        sz = ', '.join(f'{k}: {num(v)}' for k, v in sorted(sizes.items())) or '–'
        q = num(q3[c][1]) if (mo == 'sym' and c in q3) else ('–' if unc == 0 else '')
        out.append(f'| {c} {mo} | {num(fam)} | {num(unc)} | {sz} | {q} |')
out.append(f'| all | {num(totd.get("family", 0))} | {num(totd.get("open", 0))} | | |')
out.append('')
tf = totd.get('family', 0) + totd.get('open', 0)
out.append(f'So for 4 ≤ k ≤ 120 the forms of the task split into {num(totd.get("Theorem I1", 0))} covered by Theorem I1 '
           f'(k ≤ 29), {num(totd.get("family", 0))} covered by certified families ({100 * totd.get("family", 0) / tf:.1f} % '
           f'of those with k ≥ 30; 49.0 % in stage 3) and {num(totd.get("open", 0))} open, all in the sym mode; '
           f'{num(totd.get("j=k-1 (Corollary M1)", 0))} forms with j = k − 1 were set aside. Per k:')
out.append('')
out.append('| k | covered by families | open |')
out.append('|---|---|---|')
for k in (30, 40, 60, 80, 100, 120):
    if k in per:
        out.append(f'| {k} | {num(per[k].get("family", 0))} | {num(per[k].get("open", 0))} |')
out.append('')
kmins = [rows[(c, 'sym')][3] for c in ('MT', 'MU')]
out.append(f'The smallest k with an open form is {kmins[0]} in case MT and {kmins[1]} in case MU. Every open form lies '
           'in a family whose three-row part is certified (Theorem Q3), so what is open is the two-row parts Q₂ only. '
           'The checks of the script (the closed form of the presence against `tlib.norm_of` for k ≤ 30, the size ranges '
           'of Lemma Z, Lemma L, no family reported empty containing a nonempty form) ' + ('all passed.' if ok else 'did NOT all pass.'))
s = open('NOTES.md', encoding='utf-8').read()
s = s.replace('TABLE_COVERAGE', '\n'.join(out))
open('NOTES.md', 'w', encoding='utf-8').write(s)
summ = (f'**exact bookkeeping** (`coverage_st3.py`, Section 9): for 30 ≤ k ≤ 120, {num(totd.get("family", 0))} forms in '
        f'certified families ({100 * totd.get("family", 0) / tf:.1f} %) and {num(totd.get("open", 0))} open (all sym, all '
        f'with Q₃ certified); smallest open k = {min(int(x) for x in kmins)}; every anti form covered')
s = open('NOTES.md', encoding='utf-8').read()
s = s.replace('| COV_STATUS |', '| ' + summ + ' |')
open('NOTES.md', 'w', encoding='utf-8').write(s)
print('\n'.join(out))
print(summ)
