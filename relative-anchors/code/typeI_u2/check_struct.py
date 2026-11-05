"""Exact structural checks of the type I block (closed form cf.py, units of 2Y):
 (S1) the X_T-block is I - w_T w_T^T (x) E_T^{-1}  (README Remark), exactly;
 (S2) the U-variables (X_U, q) couple to X_T only through xi = X_T^T w_T: M[i, P_{t,h}] = w_t C[i, h], exactly;
 (S3) hence the Schur complement of {X_U, q, X_T} onto {X_U, q} equals the one for the 'virtual' system with a
      single T-copy of weight sqrt(s_T) (exact, at weights with s_T a rational square);
 (S4) pairwise structure: after eliminating X_T, the U-block (X_U mod kernel, q), the blocks G_tt, G_tt' of the
      Schur complement on the T-blocks (rho_t, zeta_t) do not change when the other T-weights change at fixed s_T
      (exact, Fractions).
Prints PASS."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cf import TypeI, inv_mat

ok = True
t0 = time.time()


def sub(M, I, J):
    return [[M[i][j] for j in J] for i in I]


def matmul(A, B):
    return [[sum(A[i][l] * B[l][j] for l in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def schur(M, keep, elim):
    A = sub(M, keep, keep); Bm = sub(M, keep, elim); C = sub(M, elim, elim)
    Ci = inv_mat(C, Fr)
    X = matmul(matmul(Bm, Ci), [list(r) for r in zip(*Bm)])
    return [[A[i][j] - X[i][j] for j in range(len(keep))] for i in range(len(keep))]


# (S1), (S2)
for w, T in [((Fr(1), Fr(2), Fr(3, 2), Fr(1, 2)), [0, 1]), ((Fr(2, 3), Fr(1), Fr(3, 2), Fr(1, 2), Fr(5, 4)), [1, 3]),
             ((Fr(2, 3), Fr(1), Fr(3, 2), Fr(1, 2), Fr(5, 4)), [0, 2, 4]), ((Fr(1), Fr(3), Fr(1, 2), Fr(2), Fr(1, 3)), [2])]:
    ti = TypeI(w, T); M = ti.build(); U = ti.U; s = ti.s
    sT = sum(s[t] for t in T)
    FT, ET = ti.Emat(T)
    ETi = inv_mat(ET, Fr)
    pos = {x: i for i, x in enumerate(FT)}
    e1 = max(abs(M[ti.idx[('P', t, h)]][ti.idx[('P', tp, hp)]] - ((1 if (t, h) == (tp, hp) else 0) - ti.w[t] * ti.w[tp] * ETi[pos[h]][pos[hp]]))
             for t in T for h in U for tp in T for hp in U)
    Uv = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    e2 = 0
    for i in Uv:
        for h in U:
            base = M[i][ti.idx[('P', T[0], h)]] / ti.w[T[0]]
            for t in T:
                e2 = max(e2, abs(M[i][ti.idx[('P', t, h)]] - ti.w[t] * base))
    good = e1 == 0 and e2 == 0
    ok &= good
    print(f"(S1,S2) k={len(w)} T={T}: X_T block error {e1}, coupling-through-xi error {e2}: {good}", flush=True)

# (S3): virtual copy
for wT, wU in [((Fr(3, 5), Fr(4, 5)), (Fr(1), Fr(3, 2))), ((Fr(3, 13), Fr(4, 13), Fr(12, 13)), (Fr(2), Fr(1, 2))),
               ((Fr(6, 5), Fr(8, 5)), (Fr(1), Fr(1, 3), Fr(2)))]:
    sT = sum(x * x for x in wT)
    rt = [Fr(n) for n in [3, 5, 13, 1, 2]]
    # sqrt(s_T) exactly
    num, den = sT.numerator, sT.denominator
    import math
    rn, rd = math.isqrt(num), math.isqrt(den)
    assert rn * rn == num and rd * rd == den
    wtau = Fr(rn, rd)
    m = len(wT)
    w = list(wT) + list(wU); T = list(range(m))
    ti = TypeI(w, T); M = ti.build(); U = ti.U
    keep = [ti.idx[('P', r, h)] for r in U for h in U] + [ti.idx[('Q',)]]
    elim = [ti.idx[('P', t, h)] for t in T for h in U]
    S = schur(M, keep, elim)
    w1 = [wtau] + list(wU)
    t1 = TypeI(w1, [0]); M1 = t1.build()
    keep1 = [t1.idx[('P', r, h)] for r in t1.U for h in t1.U] + [t1.idx[('Q',)]]
    elim1 = [t1.idx[('P', 0, h)] for h in t1.U]
    S1 = schur(M1, keep1, elim1)
    err = max(abs(S[i][j] - S1[i][j]) for i in range(len(keep)) for j in range(len(keep)))
    good = err == 0
    ok &= good
    print(f"(S3) m={m} u={len(wU)} s_T={sT}: U-block after eliminating X_T = virtual single copy: error {err}: {good}", flush=True)


# (S4) pairwise structure of the T-block Schur complement
def Gexact(w, T):
    ti = TypeI(w, T); M = ti.build(); U = ti.U
    r0 = U[0]
    # eliminate: X_T, X_U except the kernel coordinate P_{r0,h} (deleted), Q
    elim = [ti.idx[('P', t, h)] for t in T for h in U] + [ti.idx[('P', r, h)] for r in U if r != r0 for h in U] + [ti.idx[('Q',)]]
    keepT = {t: [ti.idx[('R', t)]] + [ti.idx[('Z', t, h)] for h in U] for t in T}
    keep = [i for t in T for i in keepT[t]]
    S = schur(M, keep, elim)
    d = len(U) + 1
    return {(T[a], T[b]): [row[b * d:(b + 1) * d] for row in S[a * d:(a + 1) * d]] for a in range(len(T)) for b in range(len(T))}


wU = [Fr(1), Fr(3, 2)]
# s_2 + s_3 = 9/4 + 1 = 13/4 = (3/5)^2 + (17/10)^2 : same s_T, different multiset of the other weights
GD = Gexact([Fr(1), Fr(2), Fr(3, 2), Fr(1)] + wU, [0, 1, 2, 3])
GE = Gexact([Fr(1), Fr(2), Fr(3, 5), Fr(17, 10)] + wU, [0, 1, 2, 3])
same01 = GD[(0, 1)] == GE[(0, 1)] and GD[(0, 0)] == GE[(0, 0)] and GD[(1, 1)] == GE[(1, 1)] and GD[(1, 0)] == GE[(1, 0)]
diff2 = GD[(2, 2)] != GE[(2, 2)]
ok &= same01 and diff2
print(f"(S4) u=2 m=4: blocks (0,0),(0,1),(1,1) unchanged when (w_2,w_3) = (3/2,1) -> (3/5,17/10) at fixed s_T: {same01}; block (2,2) changes: {diff2}", flush=True)
wU3 = [Fr(1), Fr(3, 2), Fr(1, 2)]
GF = Gexact([Fr(1), Fr(2), Fr(3, 2), Fr(1)] + wU3, [0, 1, 2, 3])
GG = Gexact([Fr(1), Fr(2), Fr(3, 5), Fr(17, 10)] + wU3, [0, 1, 2, 3])
same01b = GF[(0, 1)] == GG[(0, 1)] and GF[(0, 0)] == GG[(0, 0)] and GF[(1, 1)] == GG[(1, 1)]
diff2b = GF[(3, 3)] != GG[(3, 3)]
ok &= same01b and diff2b
print(f"(S4) u=3 m=4: blocks (0,0),(0,1),(1,1) unchanged at fixed s_T: {same01b}; block (3,3) changes: {diff2b}", flush=True)
print(f"time {time.time() - t0:.1f}s")
print("PASS" if ok else "FAIL")
