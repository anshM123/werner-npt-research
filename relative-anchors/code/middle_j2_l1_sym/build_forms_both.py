"""Collect the forms saved by cert.py in the regime (14+, 14+) (logs/both_<family>.pkl) into logs/forms_14+_14+.pkl, the
format read by check_forms.py.  Fails if a family is missing or was saved for another regime."""
import pickle, r2v
forms = {}
for F in r2v.families():
    d = pickle.load(open(f'logs/both_{F.name}.pkl', 'rb'))
    assert (d['fam'], d['s0'], d['s1']) == (F.name, '14+', '14+'), d['fam']
    forms[F.name] = d
pickle.dump(forms, open('logs/forms_14+_14+.pkl', 'wb'))
print(f'{len(forms)} families collected')
print('PASS')
