import sympy as s
h,t,r1,r2=s.symbols('h t r1 r2')
sourcejet=-1/r1+(h-1)/r2
Ceff=-1/r1+(t-1)/r2
assert s.cancel(sourcejet-Ceff-(h-t)/r2)==0
print('PASS full source0 spectator cancellation leaves (h-t)/r2')
print('PASS total 1')
