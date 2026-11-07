"""Turn a floating-point SDP solution (block coordinates, sol_*.npz) into a rational certificate file.

The certificate consists of the diagram-basis coefficient vectors of Q_A and Q_B (each 15^n numbers),
rounded to integers over a common denominator L.  Z is NOT stored: the verifier recomputes it exactly.

usage: python make_cert.py sol_file.npz n a_num a_den t_num t_den log2L out.txt
"""
import sys
import numpy as np
import blocks as bl
import brauer as br

sol, n = sys.argv[1], int(sys.argv[2])
anum, aden, tnum, tden = map(int, sys.argv[3:7])
log2L = int(sys.argv[7])
out = sys.argv[8]
pc = bl.per_copy_float()
d = np.load(sol)
yA, yB = d['yA'], d['yB']


def to_diagram(y):
    T = y.reshape([br.NDIAG] * n)
    for s in range(n):
        T = np.moveaxis(np.tensordot(pc['Phinv'], T, axes=([1], [s])), 0, s)
    return T.reshape(-1)


xA = to_diagram(yA)
xB = to_diagram(yB)
L = 2 ** log2L
XA = np.rint(xA * L).astype(object)
XB = np.rint(xB * L).astype(object)
print('max |xA|', np.abs(xA).max(), 'max |xB|', np.abs(xB).max())
print('rounding error (max) xA', np.abs(np.array(XA, dtype=float) / L - xA).max(),
      'xB', np.abs(np.array(XB, dtype=float) / L - xB).max())
with open(out, 'w') as f:
    f.write(f'n {n}\n')
    f.write(f'a {anum} {aden}\n')
    f.write(f't {tnum} {tden}\n')
    f.write(f'L {L}\n')
    for k, m in enumerate(br.DIAGRAMS):
        f.write(f'D{k} ' + ' '.join(f'{p} {q}' for p, q in m) + '\n')
    f.write('QA ' + ' '.join(str(int(v)) for v in XA) + '\n')
    f.write('QB ' + ' '.join(str(int(v)) for v in XB) + '\n')
print('wrote', out)
