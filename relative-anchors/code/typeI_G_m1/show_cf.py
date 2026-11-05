import pickle, sys, sympy as sp
u = int(sys.argv[1])
D = pickle.load(open(f'closed_t_u{u}.pkl', 'rb'))
names = D['names']
syms = sp.symbols(names)
def topoly(d):
    return sum(sp.Integer(int(c.p)) / sp.Integer(int(c.q)) * sp.prod([x**e for x, e in zip(syms, m)]) for m, c in d.items()) if d else sp.Integer(0)
G = D['G']
n = len(G)
for i in range(n):
    for j in range(i, n):
        num = sp.factor(topoly(G[i][j][0])); den = sp.factor(topoly(G[i][j][1]))
        print(f'G[{i},{j}] = ({num}) / ({den})')
        print()
