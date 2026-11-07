"""Three-copy thickness certificate, step 3: write cert/tc3_certificate.json.xz (exact, compact).

Contents (all exact; rationals as "p/q" strings, integers as ints):
  diagrams k = 2, 3, 4: canonical matchings, in the order used for every coefficient below;
  for each algebra k: integer m = 0 spin bases V_j (columns) and a right inverse Psi of the Gram-block map
     Phi[(j, r, c), a] = (V_j^T D_a V_j)[r, c]  (block coordinates in the order (j, r, c), j and r, c increasing);
  for each term (Z0: algebras 3,3,3; Z1: 4,3,3; L: 3,2,2) and spin sector beta: an integer matrix R and a rational
     symmetric Y; the operator's Gram block is R Y R^T and its diagram coefficients are (x)_s Psi_s applied to the block
     coordinates;  relations: 14 integer vectors over the 105 four-strand diagrams and the rational combination used.
The verifier (verify_tc3.py) decodes the operators and checks everything independently."""
import json
import lzma
import pickle
import numpy as np
from exact_alg import build

S = pickle.load(open('tc3_setup.pkl', 'rb'))
C = pickle.load(open('tc3_cert_blocks.pkl', 'rb'))
assert C['ok']
A = {2: build(2), 3: build(3), 4: build(4)}
alg_of = {'A2': 2, 'A3': 3, 'A4': 4}
out = {'claim': 'f(x,y) >= (1/64)(1 + (98/27) sum_s tau_s(y)) for real unit x, y on three copies, P_{3/4}',
       'strands': 'per copy (x, y, y\', w) = (0, 1, 2, 3); top points 0..k-1, bottom k..2k-1',
       'diagrams': {str(k): [[list(p) for p in m] for m in A[k]['diags']] for k in A},
       'algebras': {}, 'terms': {}, 'blocks': []}
for k in A:
    out['algebras'][str(k)] = {
        'V': {str(j): A[k]['V'][j].tolist() for j in A[k]['spins']},
        'coords': [list(c) for c in A[k]['coords']],
        'Psi': [[f'{x.numerator}/{x.denominator}' for x in row] for row in A[k]['Psi']],
    }
pub = {'Z0': 'Z0', 'Z1': 'Z1', 'L1': 'L'}
for nm, names in S['terms'].items():
    out['terms'][pub[nm]] = [alg_of[x] for x in names]
for nm in S['terms']:
    for B in S['blocks'][nm]:
        key = f'{nm}_{"".join(map(str, B["beta"]))}'
        if key not in C['Y']:
            continue
        out['blocks'].append({'term': pub[nm], 'beta': list(B['beta']), 'R': [[int(x) for x in row] for row in B['R']],
                              'Y': [[f'{p}/{q}' for (p, q) in row] for row in C['Y'][key]]})
rel_comb = []
for col, (p, q) in zip(C['rel_cols'], C['rho']):
    r_, rem = divmod(col, 225)
    a1, a2 = divmod(rem, 15)
    rel_comb.append([r_, a1, a2, f'{p}/{q}'])
out['relations'] = A[4]['relations']
out['relation_combination'] = 'sum of c * (relation r at copy 0) (x) (3-strand diagram a1, a2 at copies 1, 2, identity on w)'
out['relation_terms'] = rel_comb
with lzma.open('cert/tc3_certificate.json.xz', 'wt', encoding='utf-8') as f:
    json.dump(out, f)
print('blocks', len(out['blocks']), 'relation terms', len(rel_comb))
