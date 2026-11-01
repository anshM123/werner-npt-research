"""Negative controls for verify4m.py: writes two corrupted copies of cert4m.pkl.
  cert4m_neg_y.pkl   one multiplier changed: y_j -> y_j + 1/1024 for the first j with e_j != 0
  cert4m_neg_S.pkl   one dual diagonal entry negated (first diagonal entry of the largest dual block)
Each must be rejected by  python verify4m.py --cert <file> p_ent.
usage: python negctl4m.py"""
import pickle
from gmpy2 import mpq

def main():
    lin = pickle.load(open('lin4m.pkl', 'rb')); cert = pickle.load(open('cert4m.pkl', 'rb'))
    j = next(i for i, v in enumerate(lin['e']) if v != 0)
    y = list(cert['y']); yj = mpq(*y[j]) + mpq(1, 1024); y[j] = (int(yj.numerator), int(yj.denominator))
    pickle.dump(dict(y=y, S=cert['S']), open('cert4m_neg_y.pkl', 'wb'))
    k = max(cert['S'], key=lambda kk: len(cert['S'][kk]))
    S = {kk: [list(row) for row in M] for kk, M in cert['S'].items()}
    a, d = S[k][0][0]; S[k][0][0] = (-a, d)
    pickle.dump(dict(y=cert['y'], S=S), open('cert4m_neg_S.pkl', 'wb'))
    print(f"cert4m_neg_y.pkl: y[{j}] (e_j = {lin['e'][j]}) changed by 1/1024; cert4m_neg_S.pkl: S[{k}][0][0] = "
          f"{float(mpq(a, d)):.3e} negated (block size {len(S[k])})", flush=True)

if __name__ == "__main__":
    main()
