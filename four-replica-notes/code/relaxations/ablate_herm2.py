import sys, numpy as np, time
from alg4 import build_algebras
from herm2 import Herm2
n = 3
algs = build_algebras()
rel = Herm2(n, algs, ['SWAP', 'ORTH', 'MARG', 'PT', 'COV', 'LE2'], verbose=False)
allb = list(rel.blocks)
def onecut(name):  # PT with a single transposed replica: name like PTuuuv(3,)
    return name.startswith('PT') and name.count(',') == 1 and name.endswith(',)')
variants = {
    'all': lambda b: True,
    'no COV': lambda b: not b[0].startswith('COV'),
    'no LE2': lambda b: not b[0].startswith('LE2'),
    'COV only (+X>=0, equalities)': lambda b: b[0].startswith('COV'),
    'COV + 1-replica PT': lambda b: b[0].startswith('COV') or onecut(b[0]),
}
for vname, keep in variants.items():
    rel.blocks = [b for b in allb if keep(b)]
    for (mu, nu, s, ref) in [(0.5, 0.5, -1, 0.0), (0.8, 0.2, -1, 0.045), (0.6, 0.4, 1, 0.005)]:
        sol = rel.solve(mu, nu, s)
        print(f"[{vname}] C={mu}u{'+' if s > 0 else '-'}{nu}v: status {sol.status}, q relaxation {sol.obj_val:.6e} (true min {ref})", flush=True)
