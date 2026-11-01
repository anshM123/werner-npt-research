"""Stress adversary for BM: minimise the PPT margin of the best one-atom decomposition.
Target X = normalised(rho/Tr rho + w G^(C, D)) (mode 'nf') or normalised(G^(A1, B1) + w G^(A2, B2)) (mode 'gen'),
G^ = G/Tr G, w = exp(tau).  Objective nu(X) = mu*(X) / lambda_min(X), with
mu*(X) = max_{E, F, lam} min(lambda_min(X - lam G^(E, F)), lambda_min((X - lam G^(E, F))^Gamma))  (see ppt_margin.py).
nu < 0 robustly would mean that no single atom leaves even a PPT remainder (SF fails).  The inner maximum is tracked by
warm-started L-BFGS on a soft-min, with multi-start refreshes; the outer step uses the envelope gradient.
At the end the worst targets get: mu* from many starts, and the full one-atom fit (fit_robust).
usage: python adv_mu.py mode seed outer [lr] [cond]"""
import sys, time, numpy as np, torch
from bm_core import RHO_T, RHO, G, ptmin, rand_c, D, N, pt
from bm_fit import fit_robust, norm_t, candidates
from ppt_margin import maximize_mu, atom_hat_t, pt_t, softmin, mu_hard
torch.set_num_threads(1)

mode, seed, OUTER = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
LR = float(sys.argv[4]) if len(sys.argv) > 4 else 0.02
COND = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
BETA = 3000.0
g = np.random.default_rng(seed); torch.manual_seed(seed)


def target_t(ps, tau):
    w = torch.exp(tau)
    if mode == 'nf':
        X = RHO_T / torch.trace(RHO_T).real + w * atom_hat_t(ps[0], ps[1])
    else:
        X = atom_hat_t(ps[0], ps[1]) + w * atom_hat_t(ps[2], ps[3])
    X = (X + X.conj().T) / 2
    return X / torch.trace(X).real


nmat = 2 if mode == 'nf' else 4
ps0 = [np.linalg.qr(rand_c(g, D, D))[0] + 0.3 * rand_c(g, D, D) for _ in range(nmat)]
if COND < 1:
    U, sv, Vh = np.linalg.svd(ps0[0]); sv[-1] = COND * sv[0]; ps0[0] = U @ np.diag(sv) @ Vh
ps = [torch.tensor(p, dtype=torch.complex128).requires_grad_() for p in ps0]
tau = torch.tensor(float(g.uniform(-1, 1)), requires_grad=True)
opt = torch.optim.Adam(ps + [tau], lr=LR)


def inputs_np():
    out = [(p.detach().numpy(), q.detach().numpy()) for p, q in zip(ps[0::2], ps[1::2])]
    if mode == 'nf':
        out.append((np.eye(D), np.eye(D)))
    return out


def full_restart(X, cur=None):
    st = candidates(X) + inputs_np() + [(rand_c(g, D, D), rand_c(g, D, D)) for _ in range(3)]
    if cur is not None:
        st = [cur] + st
    return maximize_mu(X, st, beta=BETA, iters=200)


with torch.no_grad():
    X = target_t(ps, tau).numpy()
m, E, Fm, lam = full_restart(X)
print(f"stress adversary mode {mode}, seed {seed}, lr {LR}, cond {COND}: start mu* {m:+.5f}, lambda_min(X) "
      f"{np.linalg.eigvalsh(X)[0]:.4f}, PT-min {ptmin(X):+.4f}", flush=True)
Et = torch.tensor(E, dtype=torch.complex128).requires_grad_(); Ft = torch.tensor(Fm, dtype=torch.complex128).requires_grad_()
sl = torch.tensor(np.log(lam / (2 - lam)), requires_grad=True)
worst = []
t0 = time.time()
for it in range(OUTER):
    # inner: re-maximise mu over (E, F, lam) at the current target (warm start)
    with torch.no_grad():
        Xc = target_t(ps, tau)
    iopt = torch.optim.LBFGS([Et, Ft, sl], lr=1, max_iter=60, history_size=20, line_search_fn='strong_wolfe',
                             tolerance_grad=1e-12, tolerance_change=1e-14)

    def cl():
        iopt.zero_grad()
        S = Xc - 2 * torch.sigmoid(sl) * atom_hat_t(Et, Ft)
        ev = torch.cat([torch.linalg.eigvalsh(S), torch.linalg.eigvalsh(pt_t(S))])
        loss = -softmin(ev, BETA); loss.backward(); return loss
    try:
        iopt.step(cl)
    except Exception:
        pass
    Xn = Xc.numpy()
    if it % 15 == 7:
        cur = (Et.detach().numpy().copy(), Ft.detach().numpy().copy())
        m2, E2, F2, l2 = full_restart(Xn, cur)
        with torch.no_grad():
            Et.copy_(torch.tensor(E2)); Ft.copy_(torch.tensor(F2)); sl.copy_(torch.tensor(np.log(l2 / (2 - l2))))
    lamv = float(2 * torch.sigmoid(sl))
    mh = mu_hard(Xn, Et.detach().numpy(), Ft.detach().numpy(), lamv)
    lmin = np.linalg.eigvalsh(Xn)[0]
    worst.append((mh / lmin, it, [p.detach().numpy().copy() for p in ps], float(tau)))
    worst.sort(key=lambda z: z[0]); worst = worst[:5]
    # outer: descend nu = mu / lambda_min(X) with (E, F, lam) fixed (envelope gradient of the soft-min)
    X = target_t(ps, tau)
    S = X - 2 * torch.sigmoid(sl).detach() * atom_hat_t(Et.detach(), Ft.detach())
    ev = torch.cat([torch.linalg.eigvalsh(S), torch.linalg.eigvalsh(pt_t(S))])
    nu = softmin(ev, BETA) / torch.linalg.eigvalsh(X)[0]
    opt.zero_grad(); nu.backward(); opt.step()
    with torch.no_grad():
        tau.clamp_(-4.0, 4.0)
        for a, b in zip(ps[0::2], ps[1::2]):
            na, nb = torch.linalg.norm(a), torch.linalg.norm(b); a.mul_(torch.sqrt(nb / na)); b.mul_(torch.sqrt(na / nb))
    if it % 10 == 0 or it == OUTER - 1:
        sv = [np.linalg.svd(p.detach().numpy(), compute_uv=False) for p in ps]
        print(f"step {it}: mu {mh:+.5f}, nu {mh / lmin:+.4f}, lam {lamv:.3f}, w {float(torch.exp(tau)):.3f}, PT-min {ptmin(Xn):+.4f}, "
              f"lambda_min(X) {lmin:.4f}, min filter cond {min(s_[-1] / s_[0] for s_ in sv):.2e}; worst nu {worst[0][0]:+.4f} "
              f"(step {worst[0][1]})  [{time.time() - t0:.0f}s]", flush=True)

print("worst targets: mu* from many starts, then the full one-atom fit", flush=True)
res = []
for k, (nu_w, it, pp, tw) in enumerate(worst[:3]):
    with torch.no_grad():
        Xw = target_t([torch.tensor(p) for p in pp], torch.tensor(tw)).numpy()
    ins = [(pp[0], pp[1])] + ([(pp[2], pp[3])] if mode == 'gen' else [(np.eye(D), np.eye(D))])
    st = candidates(Xw) + ins + [(rand_c(g, D, D), rand_c(g, D, D)) for _ in range(10)]
    m, E, Fm, lam = maximize_mu(Xw, st, beta=BETA, iters=400)
    out = fit_robust(Xw, P=80, restarts=3, seed=2000 + k, rounds=25, extra=ins + [(E, Fm)])
    print(f"  worst #{k} (step {it}, tracked nu {nu_w:+.4f}): mu* >= {m:+.5f} (nu {m / np.linalg.eigvalsh(Xw)[0]:+.4f}), "
          f"PT-min {ptmin(Xw):+.4f}, lambda_min {np.linalg.eigvalsh(Xw)[0]:.4f}; one-atom fit: param {out['param']:.2e}, "
          f"VP {out['vp']:.2e} (atom wt {out['weight']:.3f})", flush=True)
    res.append((m, out['vp'], pp, tw, E, Fm, lam))
np.savez(f"advmu_{mode}_s{seed}.npz", mode=mode, mu=np.array([r[0] for r in res]), fit=np.array([r[1] for r in res]),
         params=np.array([np.stack(r[2]) for r in res]), tau=np.array([r[3] for r in res]),
         E=np.array([r[4] for r in res]), F=np.array([r[5] for r in res]), lam=np.array([r[6] for r in res]))
print("saved", flush=True)
