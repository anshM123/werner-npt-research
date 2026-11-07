"""eigenvalues of the blocks of Q_A and Z for a solution (diagram coordinates), float"""
import sys, itertools, numpy as np
import brauer as br, blocks as bl
d = np.load(sys.argv[1]); n = int(sys.argv[2]); a = float(sys.argv[3]); t = float(d['t'])
pc = bl.per_copy_float(); nm = pc['named']
NB = 15
piA = br.pt_perm([0]); piB = br.pt_perm([1])
S = nm['S_BB']
muL = [br.MULT[S][x][0] for x in range(NB)]; muR = [br.MULT[x][S][0] for x in range(NB)]
def perm(X, pi):
    Y = np.zeros_like(X); Y[np.ix_(*([np.array(pi)] * n))] = X; return Y
def Pc(X):
    SX = perm(X, muL); return (X + SX + perm(X, muR) + perm(SX, muR)) / 4
c = np.zeros(NB); c[nm['I']] = a; c[nm['F_AB']] -= .25; c[nm['Om_AB']] -= .25
e = np.zeros(NB); e[nm['I']] = 1
def tp(v):
    X = np.array([1.0])
    for _ in range(n): X = np.kron(X, v)
    return X.reshape([NB] * n)
xA = d['xA'].reshape([NB] * n); xB = d['xB'].reshape([NB] * n)
QA = Pc(xA)
Z = Pc(tp(c) - t * tp(e) - perm(QA, piA) - perm(xB, piB))
U = pc['U']
def blocks_of(X):
    out = {}
    for beta in itertools.product(range(4), repeat=n):
        Ub = np.array([[1.0]])
        for j in beta: Ub = np.kron(Ub, U[j])
        # operator applied via Kronecker of compressed diagram blocks
        G = 0
        for idx in itertools.product(range(NB), repeat=n):
            if X[idx] == 0: continue
            M = np.array([[1.0]])
            for s, j in enumerate(beta): M = np.kron(M, U[j].T @ pc['Dm'][idx[s]] @ U[j])
            G = G + X[idx] * M
        out[beta] = G
    return out
for name, X in [('Q_A', QA), ('Z', Z), ('Q_B', xB)]:
    if np.abs(X).max() == 0: print(name, 'zero'); continue
    bl_ = blocks_of(X)
    small = []
    for beta, G in bl_.items():
        w = np.linalg.eigvalsh((G + G.T) / 2)
        small += [(round(v, 12), beta) for v in w if abs(v) < 1e-6]
    allw = np.concatenate([np.linalg.eigvalsh((G + G.T) / 2) for G in bl_.values()])
    print(name, 't =', t, ' #|eig|<1e-6:', len(small), ' min eig', allw.min(), ' smallest |eig| above 1e-6:', np.sort(np.abs(allw[np.abs(allw) > 1e-6]))[:3])
    print('   blocks with near-zero eigenvalues:', sorted(set(b for _, b in small)))
