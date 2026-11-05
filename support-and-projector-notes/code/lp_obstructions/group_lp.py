"""Certificate search for PSD rank-two positivity at odd n using
 (i) Wu--Zou cut inequalities z_A + sqrt(e_A e'_A) - z_{AuK} - z_{AuL} >= 0 (A,K,L nonempty partition),
 (ii) grouped lower-copy theorems: for every partition B of [n] into an odd number r < n of blocks (r = 3, 5, ...),
      Z^{(B)} >= -sqrt(Qh^{(B)}(u) Qh^{(B)}(v)), valid in every local dimension (r = 3: Wu--Zou Thm 1.3;
      r >= 5 only if the scheme is used inductively).
Identity sought: Z_n = sum_B theta_B Z^{(B)} + sum_c m_c (z_A - z_{AK} - z_{AL}) + sum_T gamma_T z_T,
 theta, m, gamma >= 0, sum theta = 1, and 1 + sum_B theta_B (Qh^{(B)} - 1) + sum_c m_c e_{A_c} <= Qh_n
 coefficientwise in the local-parity distribution a_S(w) of w (x) w."""
import itertools, sys
import numpy as np
from scipy.optimize import linprog
from cut_lp import cuts


def set_partitions(items, r):
    items = list(items)
    if r == 1:
        yield [frozenset(items)]; return
    if len(items) < r:
        return
    first, rest = items[0], items[1:]
    # first in its own block
    for p in set_partitions(rest, r - 1):
        yield [frozenset([first])] + p
    # first joins an existing block
    for p in set_partitions(rest, r):
        for i in range(len(p)):
            yield p[:i] + [p[i] | {first}] + p[i + 1:]


def zcoef_block(blocks):
    """coefficients of Z^{(B)} on z_T (T a union of blocks, T != empty, != [n]) and the a-coefficients of Qh^{(B)}-1."""
    r = len(blocks)
    zc = {}
    for k in range(1, r):
        for sub in itertools.combinations(range(r), k):
            T = frozenset().union(*[blocks[i] for i in sub])
            zc[T] = zc.get(T, 0) + (-1) ** (r - k) * 2 ** k
    return zc


def qhat_minus1_coef_on_a(groups, n, S):
    """Qh^{(B)}(w) - 1 = -sum_{T union of blocks} (-1)^{r-|T|} 2^{|T|} e_T(w), e_T = 2 sum_{S:|S cap T| odd} a_S.
    Returns the coefficient of a_S."""
    r = len(groups)
    tot = 0.0
    for k in range(1, r):
        for sub in itertools.combinations(range(r), k):
            T = frozenset().union(*[groups[i] for i in sub])
            if len(S & T) % 2 == 1:
                tot += -(-1) ** (r - k) * 2 ** k * 2
    return tot


def solve(n, block_sizes=(3,)):
    U = frozenset(range(n))
    C = cuts(n)
    parts = []
    for r in block_sizes:
        parts += [p for p in set_partitions(range(n), r)]
    subsets = [frozenset(s) for k in range(1, n) for s in itertools.combinations(range(n), k)]
    nb, nc = len(parts), len(C)
    nv = nb + nc
    Zn = {T: (-1) ** (n - len(T)) * 2 ** len(T) for T in subsets}
    Aub, bub, Aeq, beq = [], [], [], []
    # gamma_T = Zn[T] - sum_B theta_B zc_B[T] - sum_c m_c (1[A=T] - 1[AuK=T] - 1[AuL=T]) >= 0
    zcs = [zcoef_block(p) for p in parts]
    for T in subsets:
        row = np.zeros(nv)
        for bi, zc in enumerate(zcs):
            row[bi] = zc.get(T, 0)
        for ci, (A, K, L) in enumerate(C):
            row[nb + ci] = (A == T) - (A | K == T) - (A | L == T)
        Aub.append(row); bub.append(Zn[T])
    row = np.zeros(nv); row[:nb] = 1; Aeq.append(row); beq.append(1.0)
    # a-condition: for every even nonempty S:  sum_B theta_B coefB(S) + sum_c m_c 2[|S cap A| odd] <= 3^|S| - 1
    for k in range(2, n + 1, 2):
        for S in itertools.combinations(range(n), k):
            S = frozenset(S)
            row = np.zeros(nv)
            for bi, p in enumerate(parts):
                row[bi] = qhat_minus1_coef_on_a(p, n, S)
            for ci, (A, K, L) in enumerate(C):
                row[nb + ci] = 2 * (len(S & A) % 2)
            Aub.append(row); bub.append(3 ** len(S) - 1)
    res = linprog(np.zeros(nv), A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(0, None)] * nv, method='highs')
    return res, parts, C


if __name__ == '__main__':
    n = int(sys.argv[1])
    sizes = tuple(int(a) for a in sys.argv[2:]) or (3,)
    res, parts, C = solve(n, sizes)
    print(f'n={n} block counts {sizes}: partitions {len(parts)}, cuts {len(C)}, status {res.status}: {res.message}')
    if res.status == 0:
        th = res.x[:len(parts)]
        for p, t in zip(parts, th):
            if t > 1e-9:
                print('  theta', round(t, 6), [sorted(b) for b in p])
        m = res.x[len(parts):]
        for c, v in zip(C, m):
            if v > 1e-9:
                print('  cut', round(v, 6), sorted(c[0]), sorted(c[1]), sorted(c[2]))
