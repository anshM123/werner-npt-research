"""Exact comparison of two independent implementations of the two-level invariant reduction, case one (u1 = 1,
u2 >= 2): the SymPy expression-tree version two_level.py (-> two_level_one_exprtree.pkl) and the polynomial-ring version
two_level_fast.py (-> two_level_one.pkl).  They share only the closed form of the block; the structured algebra (sums
over the levels, inner products) is written twice.  (The expression-tree version of the standard part,
two_level_std.py, did not finish within 400 s and is not compared; logs/abandoned/.)"""
import pickle, flint

A = pickle.load(open('two_level_one_exprtree.pkl', 'rb'))
B = pickle.load(open('two_level_one.pkl', 'rb'))
assert A['names'] == B['names'], (A['names'], B['names'])
ctx = flint.fmpq_mpoly_ctx.get(tuple(A['names']), 'lex')
ok = True
for i in range(3):
    for j in range(3):
        an, ad = ctx.from_dict(A['Q3'][i][j][0]), ctx.from_dict(A['Q3'][i][j][1])
        bn, bd = ctx.from_dict(B['Q3'][i][j][0]), ctx.from_dict(B['Q3'][i][j][1])
        ok &= (an * bd - bn * ad) == 0
print('invariant part, case one (u1 = 1): the two implementations agree exactly:', ok)
print('PASS' if ok else 'FAIL')
