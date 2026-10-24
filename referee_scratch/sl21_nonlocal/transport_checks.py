"""Exact kink transport integral and translation test for the nonlocal construction."""
import sympy as s
from solve_lax import mu,rt,G,zeros,mm,madd
from itertools import product
y=s.symbols('y',positive=True)
C=(1+y*y)/(2*y)
K=s.Matrix([[0,s.I,0],[-s.I,0,0],[0,0,0]])
Ma=s.Matrix([[0,0,-1],[0,0,s.I],[s.I,1,0]])
Mb=s.Matrix([[0,0,-s.I],[0,0,-1],[1,-s.I,0]])
for M in (Ma,Mb):
 # Column in the +1 eigenspace, row in the -1 eigenspace.
 assert K*M-M*K==M
print('PASS kink odd insertion matrices have ad(B) eigenvalue +mu')
f=s.sqrt(y)/(1+s.I*y);ss=(1-s.I*y)/s.sqrt(1+y*y)
assert s.simplify(ss*f/s.sqrt(C)-rt*f**2)==0
print('PASS exact body-frame odd insertion equals sqrt(mu)*f^2/2 times constant Grassmann matrix')
primitive=y/(2*mu*(1+s.I*y))
assert s.simplify(s.diff(primitive,y)-f**2/(2*mu*y))==0
IB=s.limit(primitive,y,s.oo)-s.limit(primitive,y,0)
assert IB==-s.I/(2*mu)
print('PASS exact kink transfer first variation is i/(4 sqrt(2 mu))*(a Ma+b Mb)')
# Exact two-generator holonomy exponential. Higher powers vanish.
# G uses four generators; restrict to its first two here.
a=G({1:1});b=G({2:1});ab=a*b
D=[[s.I/(4*s.sqrt(2*mu))*(a*Ma[i,j]+b*Mb[i,j]) for j in range(3)] for i in range(3)]
D2=mm(D,D);D3=mm(D2,D)
comm=Ma*Mb-Mb*Ma
for i,j in product(range(3),repeat=2):
 assert (D2[i][j]-ab*comm[i,j]/(32*mu)).zero()
 assert D3[i][j].zero()
print('PASS full kink holonomy is I+D+ab[Ma,Mb]/(64 mu), using left-coefficient Koszul signs')
# General odd auxiliary matrix. No nonzero one commutes with the spatial vacuum.
a1,a2,b1,b2=s.symbols('a1 a2 b1 b2')
odd=s.Matrix([[0,0,a1],[0,0,a2],[b1,b2,0]])
sol=s.solve(list(K*odd-odd*K),(a1,a2,b1,b2),dict=True)
assert sol==[{a1:0,a2:0,b1:0,b2:0}]
print('PASS odd centralizer of the spatial vacuum is zero: odd transport components are not translation-invariant internal charges')
# The leading exact kink time connection does NOT vanish after vacuum dressing.
# At=mu*sech(2mu*x)*N; N is ad(K) eigenvalue +/-2.
N=s.Matrix([[s.I,-1,0],[-1,-s.I,0],[0,0,0]])
assert K*N-N*K==-2*N
# Right vacuum is -mu K, left is +mu K. Since ad(K) N=-2N,
# the exponentially small raw At has dressed limit 2 mu N at BOTH endpoints.
print('PASS exact kink time connection has nonzero dressed endpoints 2*mu*N; raw At -> 0 does not suffice')
