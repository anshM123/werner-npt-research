"""Adversary against the normal-form ('vdm') merging rule.
Rule: bring X to its local normal form Y = (a (x) b) X (a (x) b)^* (both marginals proportional to I; Sinkhorn-type
iteration), take the eigenvector y of Y^Gamma for the smallest eigenvalue, W = polar factor of mat(y), and use the atom
G(E, F) with E = a^-1 W, F = b^-1; the weight lam maximises the PPT margin of S = X - lam G^(E, F).
The adversary minimises mu_rule(X) = max_lam min(lambda_min(S), lambda_min(S^Gamma)) over targets
X = normalised(rho/Tr rho + w G^(C, D)) ('nf') or normalised(G^(A1, B1) + w G^(A2, B2)) ('gen'), differentiating through
the whole rule (unrolled normal-form iteration, eigenvector, polar factor); lam is held at its optimum (envelope).
Every 10 steps: the best margin over all atoms mu* (maximize_mu) and the product fit of the rule's remainder.
usage: python adv_vdm.py mode seed outer [lr]"""
import sys, time, numpy as np, torch
from bm_core import RHO_T, RHO, G, ptmin, rand_c, D, N, pt, vecr, SepDict, nnls_fit, best_products, unvecr
from bm_fit import norm_t, candidates
from ppt_margin import maximize_mu, atom_hat_t, pt_t, softmin
torch.set_num_threads(1)
mode, seed, OUTER = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
LR = float(sys.argv[4]) if len(sys.argv) > 4 else 0.02
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


def msqrt_inv(R, iters=40):
    """R^{-1/2} by the coupled Newton-Schulz iteration (matrix products only, so autograd is safe)"""
    c = torch.trace(R).real; A = R / c; Y = A; Z = torch.eye(D, dtype=A.dtype)
    for _ in range(iters):
        T = 0.5 * (3 * torch.eye(D, dtype=A.dtype) - Z @ Y); Y = Y @ T; Z = T @ Z
    return Z / torch.sqrt(c)


def polar_t(M, iters=40):
    Q = M / torch.linalg.norm(M)
    for _ in range(iters):
        Q = 0.5 * (Q + torch.linalg.inv(Q).conj().T)
    return Q


def rule_atom(X, iters=60):
    a = torch.eye(D, dtype=torch.complex128); b = torch.eye(D, dtype=torch.complex128); Y = X
    for _ in range(iters):
        Y4 = Y.reshape(D, D, D, D)
        rA = torch.einsum('ijkj->ik', Y4); rB = torch.einsum('ijil->jl', Y4)
        sa, sb = msqrt_inv(rA), msqrt_inv(rB)
        K = torch.kron(sa, sb); Y = K @ Y @ K.conj().T; Y = Y / torch.trace(Y).real
        a = sa @ a; b = sb @ b
    YG = pt_t(Y); wn, Un = np.linalg.eigh(YG.detach().numpy())
    sig = wn[0] - 1e-3 * (wn[1] - wn[0]); y = torch.tensor(Un[:, 0])
    for _ in range(3):                                    # inverse iteration (differentiable through the solves)
        y = torch.linalg.solve(YG - sig * torch.eye(N, dtype=YG.dtype), y); y = y / torch.linalg.norm(y)
    k = int(torch.argmax(torch.abs(y.detach()))); y = y * (y[k].conj() / torch.abs(y[k]))
    Wp = polar_t(y.reshape(D, D))
    return (torch.linalg.inv(a) @ Wp).contiguous(), torch.linalg.inv(b).contiguous()


def best_lam(Xn, Gh):
    lams = np.linspace(0, 1.2, 241)
    vals = [min(np.linalg.eigvalsh(Xn - l * Gh)[0], np.linalg.eigvalsh(pt(Xn - l * Gh))[0]) for l in lams]
    i = int(np.argmax(vals)); lo, hi = lams[max(i - 1, 0)], lams[min(i + 1, len(lams) - 1)]
    for _ in range(40):
        m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        f1 = min(np.linalg.eigvalsh(Xn - m1 * Gh)[0], np.linalg.eigvalsh(pt(Xn - m1 * Gh))[0])
        f2 = min(np.linalg.eigvalsh(Xn - m2 * Gh)[0], np.linalg.eigvalsh(pt(Xn - m2 * Gh))[0])
        if f1 < f2: lo = m1
        else: hi = m2
    l = (lo + hi) / 2
    return l, min(np.linalg.eigvalsh(Xn - l * Gh)[0], np.linalg.eigvalsh(pt(Xn - l * Gh))[0])


def sep_fit(S, rounds=40):
    nS = np.linalg.norm(vecr(S)); dic = SepDict(g, 1500); x = vecr(S); rr = 1.0
    for rd in range(rounds):
        coef, _ = nnls_fit(x, dic.V); r = x - dic.V @ coef; rr = np.linalg.norm(r) / nS
        if rr < 1e-13: break
        new = best_products(unvecr(r), g, nstart=24, keep=12)
        xs = [p for v, p, q in new if v > 0]; ys = [q for v, p, q in new if v > 0]
        if not xs: break
        dic.add_many(xs, ys)
    return rr


nmat = 2 if mode == 'nf' else 4
ps = [torch.tensor(np.linalg.qr(rand_c(g, D, D))[0] + 0.35 * rand_c(g, D, D), dtype=torch.complex128).requires_grad_() for _ in range(nmat)]
tau = torch.tensor(float(g.uniform(-1, 1)), requires_grad=True)
opt = torch.optim.Adam(ps + [tau], lr=LR)
worst = []; t0 = time.time()
for it in range(OUTER):
    X = target_t(ps, tau)
    E, Fm = rule_atom(X)
    Gh = atom_hat_t(E, Fm)
    Xn = X.detach().numpy(); Ghn = Gh.detach().numpy()
    lam, mu = best_lam(Xn, Ghn)
    S = X - lam * Gh
    ev = torch.cat([torch.linalg.eigvalsh(S), torch.linalg.eigvalsh(pt_t(S))])
    obj = softmin(ev, BETA)
    worst.append((mu, it, [p.detach().numpy().copy() for p in ps], float(tau), lam))
    worst.sort(key=lambda z: z[0]); worst = worst[:4]
    if it % 10 == 0 or it == OUTER - 1:
        sv = [np.linalg.svd(p.detach().numpy(), compute_uv=False) for p in ps]
        ms, _, _, _ = maximize_mu(Xn, [(E.detach().numpy(), Fm.detach().numpy())] + candidates(Xn), iters=200)
        sf = sep_fit(Xn - lam * Ghn) if mu >= 0 else np.nan
        print(f"step {it}: rule margin {mu:+.6f} (lam {lam:.3f}), best margin mu* >= {ms:+.6f}, rule remainder product-fit {sf:.1e}, "
              f"PT-min {ptmin(Xn):+.4f}, w {float(torch.exp(tau)):.3f}, min filter cond {min(s_[-1] / s_[0] for s_ in sv):.2e}; "
              f"worst rule margin {worst[0][0]:+.6f}  [{time.time() - t0:.0f}s]", flush=True)
    opt.zero_grad(); obj.backward(); opt.step()
    with torch.no_grad():
        tau.clamp_(-4.0, 4.0)
        for p_, q_ in zip(ps[0::2], ps[1::2]):
            na, nb = torch.linalg.norm(p_), torch.linalg.norm(q_); p_.mul_(torch.sqrt(nb / na)); q_.mul_(torch.sqrt(na / nb))
print("worst targets for the rule:", flush=True)
for k, (mu, it, pp, tw, lam) in enumerate(worst[:2]):
    with torch.no_grad():
        Xw = target_t([torch.tensor(p) for p in pp], torch.tensor(tw)); E, Fm = rule_atom(Xw)
    Xn = Xw.numpy(); Ghn = atom_hat_t(E, Fm).numpy()
    ms, _, _, _ = maximize_mu(Xn, [(E.numpy(), Fm.numpy())] + candidates(Xn) + [(rand_c(g, D, D), rand_c(g, D, D)) for _ in range(6)], iters=300)
    sf = sep_fit(Xn - lam * Ghn) if mu >= 0 else np.nan
    print(f"  #{k} (step {it}): rule margin {mu:+.6f}, mu* >= {ms:+.6f}, rule remainder product-fit {sf:.1e}, PT-min {ptmin(Xn):+.4f}", flush=True)
np.savez(f"advvdm_{mode}_s{seed}.npz", mode=mode, mu=np.array([z[0] for z in worst]), params=np.array([np.stack(z[2]) for z in worst]),
         tau=np.array([z[3] for z in worst]))
