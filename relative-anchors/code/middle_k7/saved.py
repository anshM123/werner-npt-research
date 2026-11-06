"""Loads what a certificate run of cert5.py saved (logs/cert_*.pkl): the block over Q(q) and its leading minors, as
coefficient lists of numerators and denominators; evaluation at rational q."""
import pickle
from fractions import Fraction as Fr


def load(path):
    with open(path, 'rb') as fh:
        return pickle.load(fh)


def ev_poly(cs, x):
    v = Fr(0)
    for c in reversed(cs):
        v = v * x + Fr(c)
    return v


def ev_rf(nd, x):
    n, d = nd
    dv = ev_poly(d, x)
    if dv == 0:
        raise ZeroDivisionError
    return ev_poly(n, x) / dv


def block_at(data, x):
    return [[ev_rf(e, x) for e in row] for row in data['block']]


def minors_at(data, x):
    return [ev_rf(e, x) for e in data['minors']]


def cert_path(fam, a, b, k=5, prefix='logs/cert'):
    return f'{prefix}_k{k}_{fam}_{a}_{b}.pkl'
