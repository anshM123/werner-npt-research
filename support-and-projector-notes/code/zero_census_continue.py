"""Continue one start of zero_census.py (same n, seed and run index) for more L-BFGS rounds, printing the value, the
joint local ranks and the best cut after each round.  usage: python zero_census_continue.py n seed run_index rounds"""
import sys, itertools, importlib.util
import numpy as np
import torch

n, seed, idx, rounds = (int(a) for a in sys.argv[1:5])
sys.argv = ['zero_census.py', str(n), '0', str(seed)]          # load the census helpers without running any start
spec = importlib.util.spec_from_file_location('zero_census', 'zero_census.py')
zc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(zc)
rng = np.random.default_rng(seed)
for r in range(idx + 1):                                          # same random draws as the census
    params = [torch.tensor(rng.standard_normal((3 ** n, 2)), requires_grad=True) for _ in range(4)]
opt = torch.optim.LBFGS(params, lr=1, max_iter=5000, tolerance_grad=1e-15, tolerance_change=1e-18, history_size=80,
                        line_search_fn='strong_wolfe')


def closure():
    opt.zero_grad()
    v, _ = zc.value(params)
    v.backward()
    return v


for k in range(rounds):
    opt.step(closure)
    v, C = zc.value(params)
    C = C.detach().numpy()
    jr, _ = zc.joint_ranks(C)
    best = min(((I, zc.cut_ratio(C, I)) for kk in range(1, n) for I in itertools.combinations(range(n), kk)),
               key=lambda z: z[1])
    print(f"round {k}: q/|C|^2 = {v.item():.2e}  joint local ranks {jr}  best cut {best[0]} ratio {best[1]:.1e}",
          flush=True)
