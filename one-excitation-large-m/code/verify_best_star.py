"""Re-evaluate the best star point found at level m with the literal README definitions and q_n on (C^3)^{(x)(m+1)}."""
import glob, json, torch, math, sys
from core import *
from stage45 import fields
m = int(sys.argv[1]); N = 2 ** m
best = None
for fn in glob.glob(f'results/s45_star*_{m}_*.json'):
    d = json.load(open(fn))
    if fn.endswith('_s1.json') or fn.endswith('_s2.json'):
        continue
    if best is None or d['sup_found'] > best[0]:
        best = (d['sup_found'], fn, d['kind'])
p = torch.load(best[1].replace('.json', '.pt'))
RA, RB, LA, LB = fields(p, N, best[2].endswith('1'))
def assemble(RA, RB):
    R = torch.zeros(3, N + 1, 2, dtype=CD); R[:, :N, :] = RA.reshape(3, N, 2); R[:, N, :] = RB; return R.reshape(3 * (N + 1), 2)
R = assemble(RA, RB); L = assemble(LA, LB)
X, P, Q, Y = blocks(R, L, N)
Sig = float(Sigma_full_literal(X, m)); Gam = float(Gamma_X_literal(X, m)); y = float(q1(Y))
pl = sum(float(q1(P.reshape(3, N, 3)[:, a, :])) for a in range(N)); ql = sum(float(q1(Q.reshape(3, 3, N)[:, :, a])) for a in range(N))
ratio = math.sqrt(y) * (math.sqrt(Gam) - math.sqrt(Sig)) / math.sqrt(pl * ql)
C = R @ L.conj().T
qn = float(qn_big(C, m)) if m <= 4 else float('nan')
msg = (f"m={m} best file {best[1]} (kind {best[2]}): literal ratio={ratio:.8f} (search {best[0]:.8f}); Sigma={Sig:.6g} Gamma={Gam:.6g} "
       f"q1Y={y:.6g} q~P={pl:.4e} q~Q={ql:.4e}; rank C={int(torch.linalg.matrix_rank(C))}; Y svals={[float(x) for x in torch.linalg.svdvals(Y)[:2]]}; "
       f"q_n(C)/||C||^2={qn / float(fro2(C)):.4e}; needed 2^m-2={2**m-2}")
print(msg)
open(f'verify_best_star_m{m}.log', 'w').write(msg + '\n')
