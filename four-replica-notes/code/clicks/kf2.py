"""Claim (KF2): for Z = sum_i I_i (x) X_i on (C^d)^(x)n, with every X_i an operator on the other n-1 sites that is
fully traceless (all one-site partial traces vanish), s_1(Z)^2 + s_2(Z)^2 <= (2/3)... in general (2/d)? ||Z||_F^2.
Adversarial maximisation of (s1^2 + s2^2)/||Z||_F^2 (smooth: maximise ||P Z Q||^2 over rank-2 isometries).
usage: python kf2.py n d starts"""
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
n, d, NS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); D = d ** n
def fully_traceless(X, m):
    """project an operator on m sites onto the fully traceless sector (apply id - Tr(.)I/d at every site)"""
    T = X.reshape([d] * (2 * m))
    for j in range(m):
        tr = torch.diagonal(T, dim1=j, dim2=j + m).sum(-1)          # remove site j pair
        # rebuild tr (x) I/d on site j
        tr_e = tr.unsqueeze(-1).unsqueeze(-1) * torch.eye(d, dtype=X.dtype) / d   # indices (... , j, j')
        # move the two new axes into positions j and j+m
        perm_src = list(range(tr.ndim)) + [tr.ndim, tr.ndim + 1]
        # tr has axes: sites except j (row) then sites except j (col); insert
        rows = [a for a in range(m) if a != j]
        order = []
        k = 0
        full_axes = []
        # build tensor with axes ordering rows(m) cols(m)
        tr_e = tr_e.permute(*range(tr.ndim), tr.ndim, tr.ndim + 1)
        # current: [rows without j][cols without j][j_row][j_col]
        cur = rows + [m + a for a in rows] + [j, m + j]
        inv = [cur.index(t) for t in range(2 * m)]
        T = T - tr_e.permute(*inv)
    return T.reshape(d ** m, d ** m)
def build_Z(Xs):
    Z = 0
    for i in range(n):
        X = fully_traceless(Xs[i], n - 1)
        # I on site i, X on the others
        Xt = X.reshape([d] * (2 * (n - 1)))
        full = torch.einsum('...,ab->...ab', Xt, torch.eye(d, dtype=X.dtype))   # axes: rows(n-1) cols(n-1) a b
        rows = [a for a in range(n) if a != i]
        cur = rows + [n + a for a in rows] + [i, n + i]
        inv = [cur.index(t) for t in range(2 * n)]
        Z = Z + full.permute(*inv).reshape(D, D)
    return Z
best = 0
for s in range(NS):
    torch.manual_seed(s)
    Xs = [torch.randn(d ** (n - 1), d ** (n - 1), dtype=torch.complex128, requires_grad=True) for _ in range(n)]
    U = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True); V = torch.randn(D, 2, dtype=torch.complex128, requires_grad=True)
    o = torch.optim.Adam(Xs + [U, V], lr=0.02)
    for it in range(2500):
        Z = build_Z(Xs)
        Qu, _ = torch.linalg.qr(U); Qv, _ = torch.linalg.qr(V)
        val = torch.sum(torch.abs(Qu.conj().T @ Z @ Qv) ** 2) / torch.sum(torch.abs(Z) ** 2)
        o.zero_grad(); (-val).backward(); o.step()
    with torch.no_grad():
        Z = build_Z(Xs); sv = torch.linalg.svdvals(Z); r = (sv[0] ** 2 + sv[1] ** 2) / torch.sum(torch.abs(Z) ** 2)
        # sanity: fully traceless check of the X parts: Tr_ij Z should vanish
    best = max(best, r.item())
    print(f"start {s}: (s1^2 + s2^2)/||Z||^2 = {r.item():.6f}", flush=True)
print(f"n={n} d={d}: max ratio {best:.6f}   (claim: <= 2/3 = {2/3:.6f})")
