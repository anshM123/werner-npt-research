"""Continue gradient ascent of the (star) ratio from the best point found at level m (lower lr, more steps)."""
import glob, json, torch, math, sys
from core import *
from stage45 import fields, ratio_terms
torch.set_num_threads(2)
m = int(sys.argv[1]); N = 2 ** m; steps = int(sys.argv[2])
best = None
for fn in glob.glob(f'results/s45_star_{m}_*.json'):
    if fn.endswith('_s1.json') or fn.endswith('_s2.json'):
        continue
    d = json.load(open(fn))
    if best is None or d['sup_found'] > best[0]:
        best = (d['sup_found'], fn)
p = torch.load(best[1].replace('.json', '.pt')).clone().requires_grad_(True)
hist = []
for lr in (3e-4, 1e-4, 3e-5):
    opt = torch.optim.Adam([p], lr=lr)
    for it in range(steps):
        val, info = ratio_terms(*fields(p, N, False), N, 'star', m)
        v = float(val.detach())
        if not math.isfinite(v):
            break
        opt.zero_grad(); (-val).backward(); opt.step()
    hist.append((lr, v, {k: float(t.detach()) for k, t in info.items()}))
msg = f"m={m} start {best[0]:.6f} -> " + "; ".join(f"lr {lr:g}: {v:.6f} (p={i['p']:.2e}, q={i['q']:.2e})" for lr, v, i in hist)
print(msg); open(f'continue_star_m{m}.log', 'w').write(msg + '\n')
torch.save(p.detach(), f'results/continue_star_m{m}.pt')
