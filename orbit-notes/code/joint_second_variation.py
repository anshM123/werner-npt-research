"""Check the exact joint second variation at squeeze-critical points (t = 1/2).

For C critical on its squeeze orbit, r = R(C), and a 1-local direction H = sum_s H_s (H_s at copy s),
C(eps) = exp(eps H) C exp(eps H) has  R(C(eps)) = r + eps^2 Hq(H)/||C||^2 + O(eps^3)  with

    Hq(H) = 2 q_r({H, C}) - 2 sum_s [ q_r(H_s C) + q_r(C H_s) ]
            + 8 t^2 sum_{s<s'} Re[ B''(Tr_ss' C, Tr_ss'(H_s H_s' C)) - B''(Tr_ss'(H_s C), Tr_ss'(H_s' C)) ],

where q_r = q - r ||.||^2, {H, C} = HC + CH, B'' is the form on the other n-2 copies (the product of the two traces
when n = 2).  Equivalently Hq = 4 sum_{k<l} Re B_r(V_k, V_l) + (trace terms) over the 2n rank-two pieces
V in {H_s C, C H_s}.  The script compares Hq with v^T Hess(R) v / 2 from automatic differentiation.
usage: python joint_second_variation.py   (run from orbit-notes/)"""
import sys, itertools
import numpy as np
import torch

files = [('data/critical_n2_R0.2.npy', 2), ('data/copywise_min_n3.npy', 3)]
sys.path.insert(0, 'code')
t = 0.5
for fn, n in files:
    sys.argv = ['x', str(n), '0.5', '0.1', '0', '0']
    if 'squeeze_hessian' in sys.modules:
        del sys.modules['squeeze_hessian']
    import squeeze_hessian as SH

    def emb(M, s):
        out = np.ones((1, 1), complex)
        for k in range(n):
            out = np.kron(out, M if k == s else np.eye(3))
        return out

    def ptr(Y, S, m):
        T = Y.reshape((3,) * (2 * m))
        for s in sorted(S, reverse=True):
            T = np.trace(T, axis1=s, axis2=T.ndim // 2 + s)
        k = m - len(S)
        return T.reshape(3 ** k, 3 ** k)

    def Bpp(Y, Z, m):
        return np.conj(Y.reshape(-1)[0]) * Z.reshape(-1)[0] if m == 0 else SH.b_np(Y, Z, m)

    C = np.load(fn); C = C / np.linalg.norm(C); r = SH.q_np(C, n)
    qr_ = lambda Y: SH.q_np(Y, n) - r * np.vdot(Y, Y).real
    Ct = torch.tensor(C)
    Hs = torch.autograd.functional.hessian(lambda hh: SH.squeeze_R(Ct, hh), torch.zeros(8 * n))
    Hs = ((Hs + Hs.T) / 2).numpy()
    rng = np.random.default_rng(1)
    print(f"{fn} (n = {n}, R = {r:+.6f})")
    for trial in range(3):
        v = rng.standard_normal(8 * n)
        Hl = [SH.herm8(torch.tensor(v[8 * s:8 * s + 8])).numpy() for s in range(n)]
        E = [emb(Hl[s], s) for s in range(n)]
        Hbig = sum(E)
        val = 2 * qr_(Hbig @ C + C @ Hbig) - 2 * sum(qr_(E[s] @ C) + qr_(C @ E[s]) for s in range(n))
        for s, s2 in itertools.combinations(range(n), 2):
            a = Bpp(ptr(C, {s, s2}, n), ptr(E[s] @ E[s2] @ C, {s, s2}, n), n - 2)
            b = Bpp(ptr(E[s] @ C, {s, s2}, n), ptr(E[s2] @ C, {s, s2}, n), n - 2)
            val += 8 * t * t * (a - b).real
        print(f"   random direction {trial}: formula {val:+.10f}   autograd {v @ Hs @ v / 2:+.10f}")
