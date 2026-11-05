"""Exact (integer) tensor-product Bernstein certificates of positivity on [0,1]^n with dyadic subdivision.

A polynomial is a dense numpy object array C (Python ints), C[i_1,...,i_n] = coefficient of x_1^i_1 ... x_n^i_n.
For a dyadic box prod_j [k_j / 2^L_j, (k_j + 1) / 2^L_j] we substitute x_j = (k_j + y_j) / 2^L_j and multiply by
2^(L_j d_j) (a positive factor), which keeps integer coefficients.  The scaled Bernstein coefficients
    beta_K = sum_{I <= K} prod_j C(d_j - i_j, k_j - i_j) q_I     ( = prod_j C(d_j, k_j) * b_K )
are integers with the signs of the true Bernstein coefficients b_K.  Since p = sum_K b_K B_K with B_K >= 0 and
sum_K B_K = 1 on the box, min_K b_K > 0 certifies p > 0 on the closed box.

Corner rule (used only where stated): if every beta_K > 0 except beta_{K0} = 0 at one vertex index K0, then p > 0 on
the box minus that vertex (p(x) >= min_{K != K0} b_K (1 - B_{K0}(x)) and B_{K0}(x) < 1 off the vertex).
"""
import numpy as np
from math import comb

_cache_T = {}
_cache_S = {}


def zeros(shape):
    A = np.empty(shape, dtype=object)
    A.fill(0)
    return A


def from_dict(dct, degs):
    A = zeros(tuple(d + 1 for d in degs))
    for mon, c in dct.items():
        A[tuple(mon)] += int(c)
    return A


def apply_axis(A, M, axis):
    B = np.tensordot(M, A, axes=([1], [axis]))
    return np.moveaxis(B, 0, axis)


def bern_matrix(d):
    if d not in _cache_T:
        T = zeros((d + 1, d + 1))
        for kk in range(d + 1):
            for i in range(kk + 1):
                T[kk, i] = comb(d - i, kk - i)
        _cache_T[d] = T
    return _cache_T[d]


def shift_matrix(d, L, k):
    key = (d, L, k)
    if key not in _cache_S:
        M = zeros((d + 1, d + 1))
        for r in range(d + 1):
            for i in range(r, d + 1):
                M[r, i] = comb(i, r) * k ** (i - r) * 2 ** (L * (d - i))
        _cache_S[key] = M
    return _cache_S[key]


def bernstein_on_box(C, box):
    """box: tuple of (L_j, k_j).  Returns integer array of scaled Bernstein coefficients."""
    A = C
    for j, (L, k) in enumerate(box):
        d = C.shape[j] - 1
        if L > 0 or k > 0:
            A = apply_axis(A, shift_matrix(d, L, k), j)
        A = apply_axis(A, bern_matrix(d), j)
    return A


def certify(C, max_boxes=200000, corner=None, log=None, split_rule='width'):
    """Prove C > 0 on [0,1]^n (closed), except possibly at the single point `corner` (a tuple in {0,1}^n), where
    C must vanish and the corner rule is applied.  Returns (ok, stats)."""
    n = C.ndim
    stack = [tuple((0, 0) for _ in range(n))]
    done = 0
    corner_boxes = 0
    maxL = [0] * n
    minval = None
    while stack:
        box = stack.pop()
        B = bernstein_on_box(C, box)
        flat = B.ravel()
        mn = min(flat)
        if mn > 0:
            done += 1
            continue
        if corner is not None:
            # is a vertex of this box the corner point?
            at = all((cj == 0 and k == 0) or (cj == 1 and k == 2 ** L - 1) for cj, (L, k) in zip(corner, box))
            if at:
                K0 = tuple(0 if cj == 0 else C.shape[j] - 1 for j, cj in enumerate(corner))
                if B[K0] == 0:
                    B2 = B.copy()
                    B2[K0] = 1
                    if min(B2.ravel()) > 0:
                        done += 1
                        corner_boxes += 1
                        continue
        if done + len(stack) > max_boxes:
            return False, dict(boxes=done, failed_box=box, reason='max_boxes')
        # split
        if split_rule == 'width':
            j = min(range(n), key=lambda jj: (box[jj][0], -C.shape[jj]))
        else:
            j = split_rule(box, B)
        L, k = box[j]
        if L > 40:
            return False, dict(boxes=done, failed_box=box, reason='depth')
        for kk in (2 * k, 2 * k + 1):
            nb = list(box)
            nb[j] = (L + 1, kk)
            stack.append(tuple(nb))
        maxL[j] = max(maxL[j], L + 1)
    return True, dict(boxes=done, corner_boxes=corner_boxes, maxL=maxL)
