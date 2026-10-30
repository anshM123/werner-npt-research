"""Budget with CORRELATED separable factors: four qutrit copies from <= m endpoint cores.

Generators (all one-copy undistillable given the m-copy endpoint theorem, by absorption of separable factors):
    E^{(x) C}  (x)  sigma_R ,     C = core sites (|C| <= m), R = the other sites,
where E = X(3) is the endpoint Werner operator and sigma_R = |alpha><alpha| (x) |beta><beta| is a product state across
Alice|Bob on the sites R (alpha, beta arbitrary vectors in (C^3)^{(x)|R|}, so CORRELATED across the sites of R).
After the Werner twirl at every site the generator is fixed by its sector weights; profiles factorise between C and R.
Separable profile on R: p_T = 2^{-|R|} sum_{U subset R} (-1)^{|U cap T|} Tr[rho^alpha_U rho^beta_U].
LP (S_4-symmetrised, 5 aggregates P_k) with column generation over (alpha, beta) by Adam, for each |C| = 0..m.
usage: corr_sep_lp.py beta m starts steps rounds [seed]
"""
import itertools
import sys
import time
from math import comb

import numpy as np
import torch
from scipy.optimize import linprog

torch.set_default_dtype(torch.float64)
beta = float(sys.argv[1]); m = int(sys.argv[2])
starts, steps, rounds = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
g = torch.Generator().manual_seed(seed)
E1 = np.array([0.4, 0.6])                      # endpoint per-site profile (s, a) at d = 3


def target(b):
    ws, wa = 6 * (1 - b), 3 * (1 + b)
    P = np.array([comb(4, k) * wa ** k * ws ** (4 - k) for k in range(5)])
    return P / P.sum()


def sep_profile(alpha, beta_, r):
    """alpha, beta_ : tensors of shape (3,)*r (complex), normalised.  Returns dict T(tuple) -> weight."""
    rhoA = torch.einsum(alpha.conj(), list(range(r)), alpha, list(range(r, 2 * r)), list(range(2 * r))) if False else None
    out = {}
    tr = {}
    for U in itertools.chain.from_iterable(itertools.combinations(range(r), k) for k in range(r + 1)):
        # Tr[rho^a_U rho^b_U] = sum over kept indices of <a|..|a><b|..|b> contraction
        rest = [i for i in range(r) if i not in U]
        A = alpha.reshape([3] * r); B = beta_.reshape([3] * r)
        # rho^a_U[x, x'] = sum_rest a[x, rest] conj(a[x', rest]);  Tr[rho_a rho_b] = sum rho_a[x,x'] rho_b[x',x]
        la = list(range(r)); la2 = [i if i in rest else r + i for i in range(r)]
        rhoa = torch.einsum(A, la, A.conj(), la2, [i for i in U] + [r + i for i in U])
        rhob = torch.einsum(B, la, B.conj(), la2, [i for i in U] + [r + i for i in U])
        k = len(U)
        if k == 0:
            tr[U] = (A.abs() ** 2).sum() * (B.abs() ** 2).sum()
        else:
            ra = rhoa.reshape(3 ** k, 3 ** k); rb = rhob.reshape(3 ** k, 3 ** k)
            tr[U] = torch.real(torch.trace(ra @ rb))
    for T in itertools.chain.from_iterable(itertools.combinations(range(r), k) for k in range(r + 1)):
        out[T] = sum((-1) ** len(set(U) & set(T)) * tr[U] for U in tr) / 2 ** r
    return out


def gen_profile5(c, sepP, r):
    """c = number of core sites (placed on sites 0..c-1), separable factor on sites c..3 with profile sepP (dict over
    subsets of range(r)).  Returns S_4-symmetrised aggregate P_k."""
    P = torch.zeros(5)
    for coreT in itertools.product(range(2), repeat=c):
        wcore = float(np.prod([E1[t] for t in coreT]))
        for T, v in sepP.items():
            P[sum(coreT) + len(T)] += wcore * v
    return P / P.sum()


tgt = target(beta)
cols, labels = [], []
# seed columns: plain budget products (O, X2, E per site) with <= m cores
site = {'O': np.array([1.0, 0.0]), 'E': E1, 'X2': np.array([0.5, 0.5])}
for combo in itertools.combinations_with_replacement(['O', 'X2', 'E'], 4):
    if combo.count('E') > m:
        continue
    P = np.zeros(5)
    for types in itertools.product(range(2), repeat=4):
        P[sum(types)] += np.prod([site[combo[i]][types[i]] for i in range(4)])
    cols.append(P / P.sum()); labels.append(('prod', combo))


def lp(cols):
    G = np.array(cols); n_ = len(cols)
    res = linprog(np.concatenate([-tgt, [1.0]]), A_ub=np.hstack([G, -np.ones((n_, 1))]), b_ub=np.zeros(n_),
                  bounds=[(-1, 1)] * 5 + [(None, None)], method='highs')
    return -res.fun, res.x[:5], res.x[5]


t0 = time.time()
gap, w, c0 = lp(cols)
print(f'beta={beta} m={m}: products-only gap {gap:+.6f}', flush=True)
for rd in range(rounds):
    best = (-1e9, None, None)
    wt = torch.tensor(w)
    for c in range(0, m + 1):
        r = 4 - c
        if r == 0:
            continue
        for s in range(starts):
            a = torch.randn((3,) * r, generator=g, dtype=torch.complex128) if False else \
                (torch.randn((3,) * r, generator=g) + 1j * torch.randn((3,) * r, generator=g))
            b = (torch.randn((3,) * r, generator=g) + 1j * torch.randn((3,) * r, generator=g))
            if s == 0:      # structured start: maximally entangled pairs / GHZ
                a = torch.zeros((3,) * r, dtype=torch.complex128); b = torch.zeros((3,) * r, dtype=torch.complex128)
                for k in range(3):
                    a[(k,) * r] = 1; b[(k,) * r] = 1
                a = a + 0.01 * torch.randn((3,) * r, generator=g); b = b + 0.01 * torch.randn((3,) * r, generator=g)
            a.requires_grad_(True); b.requires_grad_(True)
            opt = torch.optim.Adam([a, b], lr=0.05)
            for it in range(steps):
                opt.zero_grad()
                an = a / torch.sqrt((a.abs() ** 2).sum()); bn = b / torch.sqrt((b.abs() ** 2).sum())
                P = gen_profile5(c, sep_profile(an, bn, r), r)
                loss = -(wt @ P); loss.backward(); opt.step()
            with torch.no_grad():
                an = a / torch.sqrt((a.abs() ** 2).sum()); bn = b / torch.sqrt((b.abs() ** 2).sum())
                P = gen_profile5(c, sep_profile(an, bn, r), r)
                val = float(wt @ P)
            if val > best[0]:
                best = (val, P.numpy().copy(), (c, an.detach().numpy().copy(), bn.detach().numpy().copy()))
    cols.append(best[1]); labels.append(('sep', best[2][0]))
    gap, w, c0 = lp(cols)
    print(f'  round {rd}: best generator has {best[2][0]} cores, value {best[0]:+.5f}; gap {gap:+.6f} '
          f'({time.time() - t0:.0f}s)', flush=True)
    if gap <= 1e-9:
        break
print(f'RESULT beta={beta} m={m}: gap {gap:+.6f} ({"IN hull" if gap <= 1e-9 else "NOT in hull"})', flush=True)
np.save(f'corr_sep_cols_{beta}_{m}.npy', np.array(cols))
