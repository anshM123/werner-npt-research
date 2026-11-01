"""Negative controls for verifyL.py: deliberately corrupted copies of certL.pkl (each must FAIL verification).
  certL_neg_y.pkl : one multiplier changed, y+_0 -> y+_0 + 1/1024;
  certL_neg_S.pkl : one dual diagonal entry negated (the largest diagonal entry of the C (-,-) part of the largest
                    lifted block with both parts);
  certL_neg_B.pkl : one cross entry B_k[0, 0] of the same block changed by + 1/1024 (only the Z identity involves it).
usage: python negctlL.py"""
import pickle
from fractions import Fraction

def main():
    cert = pickle.load(open('certL.pkl', 'rb')); fd = pickle.load(open('facedataLr.pkl', 'rb'))
    blocks = fd['blocks']
    # (1) multiplier
    c1 = {'y': {g: list(v) for g, v in cert['y'].items()}, 'S': cert['S']}
    a, q = c1['y']['+'][0]; f = Fraction(a, q) + Fraction(1, 1024); c1['y']['+'][0] = (f.numerator, f.denominator)
    pickle.dump(c1, open('certL_neg_y.pkl', 'wb'))
    print(f"certL_neg_y.pkl: y+_0 = {float(Fraction(a, q)):.6e} -> {float(f):.6e}")
    # (2) negated diagonal entry and (3) perturbed cross entry, in the largest block with both parts
    k = max((k for k, b in enumerate(blocks) if b['rp'] and b['rm']), key=lambda k: blocks[k]['rp'] + blocks[k]['rm'])
    rp, rm = blocks[k]['rp'], blocks[k]['rm']
    M = [list(row) for row in cert['S'][k]]
    i = max(range(rp, rp + rm), key=lambda i: Fraction(*M[i][i]))
    a, q = M[i][i]; M[i][i] = (-a, q)
    S2 = dict(cert['S']); S2[k] = M
    pickle.dump({'y': cert['y'], 'S': S2}, open('certL_neg_S.pkl', 'wb'))
    print(f"certL_neg_S.pkl: block {blocks[k]['name']} (r+ {rp}, r- {rm}), diagonal entry {i} = {float(Fraction(a, q)):.6e} negated")
    M3 = [list(row) for row in cert['S'][k]]
    f = Fraction(*M3[0][rp]) + Fraction(1, 1024)
    M3[0][rp] = (f.numerator, f.denominator); M3[rp][0] = (f.numerator, f.denominator)
    S3 = dict(cert['S']); S3[k] = M3
    pickle.dump({'y': cert['y'], 'S': S3}, open('certL_neg_B.pkl', 'wb'))
    print(f"certL_neg_B.pkl: block {blocks[k]['name']}, cross entry B[0, 0] -> {float(f):.6e} (+ 1/1024)")

if __name__ == "__main__":
    main()
