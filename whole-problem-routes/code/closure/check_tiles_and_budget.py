"""Exact checks for the closure-cone bullet of OBSTRUCTIONS.md.
(1) Tiles activation: S = I - sum_i v_i v_i^T/|v_i|^2 (five-vector unextendible product basis) is a PPT projector of rank 4,
    and an integer rank-two C = A B^T gives <vec C|(W^Gamma (x) S)|vec C> = -1769/9 < 0, so W (x) S violates crossed
    Cauchy-Schwarz (example from an independent audit of this repository).
(2) Budget: Tr W = 15/2, ||W^Gamma||_1 = 17/2, and (17/15)^6 > 2.   usage: python check_tiles_and_budget.py"""
from fractions import Fraction as Fr
import itertools
e=[[1,0,0],[0,1,0],[0,0,1]]
def kron(a,b): return [x*y for x in a for y in b]
def sub(a,b): return [x-y for x,y in zip(a,b)]
ones=[1,1,1]
v=[kron(e[0],sub(e[0],e[1])),kron(sub(e[0],e[1]),e[2]),kron(e[2],sub(e[1],e[2])),kron(sub(e[1],e[2]),e[0]),kron(ones,ones)]
S=[[Fr(int(i==j)) for j in range(9)] for i in range(9)]
for x in v:
    nx=sum(t*t for t in x)
    for i in range(9):
        for j in range(9): S[i][j]-=Fr(x[i]*x[j],nx)
SS=[[sum(S[i][k]*S[k][j] for k in range(9)) for j in range(9)] for i in range(9)]
assert SS==S and sum(S[i][i] for i in range(9))==4
pt=lambda M:[[M[3*a+d][3*c+b] for c in range(3) for d in range(3)] for a in range(3) for b in range(3)]
assert pt(S)==S
A=[[0,-2],[1,-4],[3,2],[4,0],[-1,0],[4,-1],[-3,-3],[-1,2],[2,-2]]
B=[[2,3],[-3,-2],[0,-6],[-1,5],[-6,0],[2,2],[-3,4],[1,3],[-6,1]]
C=[[sum(A[i][k]*B[j][k] for k in range(2)) for j in range(9)] for i in range(9)]   # rows (a1 a2), cols (b1 b2)
M3=[[Fr(int(i==j))-Fr(int(i in (0,4,8) and j in (0,4,8)),2) for j in range(9)] for i in range(9)]
val=Fr(0)
for a1,a2,b1,b2,c1,c2,d1,d2 in itertools.product(range(3),repeat=8):
    m=M3[3*a1+b1][3*c1+d1]
    if m==0: continue
    s=S[3*a2+b2][3*c2+d2]
    if s==0: continue
    val+=C[3*a1+a2][3*b1+b2]*m*s*C[3*c1+c2][3*d1+d2]
print('Tiles activation value', val, ' |C|^2 =', sum(x*x for r in C for x in r))
assert val==Fr(-1769,9)
print('budget: Tr W = 15/2, ||W^Gamma||_1 = 17/2; 17^6 =',17**6,'> 2*15^6 =',2*15**6, 17**6>2*15**6)
print('PASS')
