"""Writes the family lists of famlist.py to lists/: CASE_MODE_NS.txt (NS symbolic sizes), CASE_anti_all.txt (every anti
family), CASE_sym_2e.txt (two symbolic sizes, not both i and tau': elimination), CASE_sym_2d.txt (i and tau' symbolic,
q0 and t fixed: the determinant route, one 'q0,t' per line).  Prints the counts."""
import os
import famlist

os.makedirs('lists', exist_ok=True)
for case in ('MT', 'MU'):
    for mode in ('sym', 'anti'):
        allf = []
        for ns in (1, 2, 3, 4):
            fs = famlist.families(case, mode, ns, ns)
            allf += fs
            with open(f'lists/{case}_{mode}_{ns}.txt', 'w') as fh:
                fh.write('\n'.join(famlist.spec(case, mode, f) for f in fs) + '\n')
            print(case, mode, ns, len(fs))
        if mode == 'anti':
            with open(f'lists/{case}_anti_all.txt', 'w') as fh:
                fh.write('\n'.join(famlist.spec(case, mode, f) for f in allf) + '\n')
            print(case, mode, 'all', len(allf))
        else:
            two = famlist.families(case, mode, 2, 2)
            e = [f for f in two if not (f[2] is None and f[3] is None)]
            d = [f for f in two if f[2] is None and f[3] is None]
            with open(f'lists/{case}_sym_2e.txt', 'w') as fh:
                fh.write('\n'.join(famlist.spec(case, mode, f) for f in e) + '\n')
            with open(f'lists/{case}_sym_2d.txt', 'w') as fh:
                fh.write('\n'.join(f'{f[0]},{f[1]}' for f in d) + '\n')
            print(case, mode, '2e', len(e), '2d', len(d))
print('PASS')
