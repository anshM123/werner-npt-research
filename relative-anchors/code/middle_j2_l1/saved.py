"""Loading the forms S certified by cert_fixed.py (files logs/cert_k{k}_{u0}_{u1}.pkl written by run_all.sh)."""
import os, pickle
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def poly_at(cs, x):
    v = Fr(0)
    for c in reversed(cs):
        v = v * x + Fr(c)
    return v


def load_saved(k, u0, u1):
    """dict irrep -> matrix of functions q -> S[i][j](q) (exact, Fractions)"""
    with open(os.path.join(HERE, 'logs', f'cert_k{k}_{u0}_{u1}.pkl'), 'rb') as fh:
        dat = pickle.load(fh)
    assert (dat['u0'], dat['u1'], dat['b']) == (u0, u1, '1')
    out = {}
    for lam, dim, rows in dat['parts']:
        out[tuple(tuple(x) for x in lam)] = [[(lambda nc, dc: (lambda x: poly_at(nc, x) / poly_at(dc, x)))(nc, dc)
                                              for (nc, dc) in row] for row in rows]
    return out
