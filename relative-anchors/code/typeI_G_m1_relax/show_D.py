"""Display the local blocks D(s) of the reduced system symbolically (which constants they involve)."""
import pickle, sympy as sp
from reduced_eval import eval_srepr
D = pickle.load(open('reduced.pkl', 'rb'))
s, sig, S = sp.symbols('s sigma S', positive=True)
IG, IGT, u = sp.symbols('IG IGT u')
for (a, b), v in D['loc'].items():
    e = eval_srepr(v)
    fs = sorted(str(x) for x in e.free_symbols)
    e2 = sp.factor(sp.cancel(sp.together(e)))
    print((a, b), 'symbols:', fs)
    print('   ', e2)
