"""Find zeros of q_4(P_u - P_v) (t = -1) by local search and save U = [u, v] (81 x 2, complex).
usage: python tm1_zeros.py seed [iters]   -> tm1z/U_<seed>.npy and a one-line report"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from local_herm import run

if __name__ == "__main__":
    seed = int(sys.argv[1]); iters = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    n = int(os.environ.get('TM1_N', '4')); t = float(os.environ.get('TM1_T', '-1'))
    os.makedirs('tm1z', exist_ok=True)
    t0 = time.time()
    best, U = run(n, t, restarts=1, iters=iters, seed=seed)
    tag = f"n{n}_t{t:+.0f}_s{seed}"
    np.save(f'tm1z/U_{tag}.npy', U)
    print(f"{tag}: best q {best:.3e}; {time.time() - t0:.0f}s", flush=True)
