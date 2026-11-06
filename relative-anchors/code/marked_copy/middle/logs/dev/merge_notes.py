"""development helper (not part of any proof): fills the placeholders of NOTES.md from NOTES_body.md, NOTES_sec6.md
and a file of values (lines KEY=VALUE, '\\n' for a newline)"""
import sys, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
vals = {}
for line in open(sys.argv[1], encoding='utf-8'):
    line = line.rstrip('\n')
    if '=' in line:
        k, v = line.split('=', 1)
        vals[k] = v.replace('\\n', '\n')
notes = open('NOTES.md', encoding='utf-8').read()
notes = notes.replace('{{BODY}}', open('NOTES_body.md', encoding='utf-8').read())
notes = notes.replace('{{SEC6}}', open('NOTES_sec6.md', encoding='utf-8').read())
for k, v in vals.items():
    notes = notes.replace('{{' + k + '}}', v)
left = sorted({p.split('}}')[0] for p in notes.split('{{')[1:]})
print('UNFILLED:', left)
open('NOTES.md', 'w', encoding='utf-8').write(notes)
print('written', len(notes))
