"""Compressed products Y(k,s,t) = (Pi x Pi)(I^k (x) A^s (x) B^t)(Pi x Pi), A = P_an, B = A^Gamma = (I - Omega Omega^T)/2."""
import numpy as np
from zbasis import one_copy_ops, transfer

def build_Y(n, d=3, ops_list=None):
    ops = one_copy_ops(d)
    if ops_list is None:
        A = ops['Pan']; B = (ops['I'] - ops['Om'])/2
        ops_list = [ops['I'], A, B]
    Y = {(0, 0, 0): np.array([[1.0]])}
    for m in range(1, n+1):
        Ynew = {}
        for e1 in range(m+1):
            for e2 in range(m+1-e1):
                e3 = m - e1 - e2
                e = (e1, e2, e3)
                for tpos in range(3):
                    if e[tpos] > 0:
                        prev = list(e); prev[tpos] -= 1
                        Ynew[e] = transfer(Y[tuple(prev)], m, ops_list[tpos], d)
                        break
        Y = Ynew
    return Y
