import sympy as s
x,y,u,v,q0,q2=s.symbols('x y u v q0 q2')
p=s.Function('psi')(x,y);A=s.Function('A')(x,y);B=s.Function('B')(x,y);H=s.Function('H')(x,y)
# Complete degree-one polarization in each pair. H is both same-chirality Wick crossmatchings.
G=p*(1+u*A+v*B+u*v*(A*B+H))
def M0(f):return q0*q0*s.diff(f,x)-q0*s.diff(f,u)
def M2(f):return q2*q2*s.diff(f,y)-q2*s.diff(f,v)
val=M0(M2(G)).subs({u:0,v:0})
rev=M2(M0(G)).subs({u:0,v:0})
assert s.expand(val-rev)==0
print('PASS unequal-source crossing diamond on complete Gaussian generator')
want=q0*q0*q2*q2*s.diff(p,x,y)-q0*q0*q2*s.diff(B*p,x)-q0*q2*q2*s.diff(A*p,y)+q0*q2*(A*B+H)*p
assert s.expand(val-want)==0
print('PASS full joint first jets and both derivative Wick crossmatchings retained')
R,d0,d2,d02=s.symbols('R delta0 delta2 delta02')
R00=R;R10=R+d0;R01=R+d2;R11=R+d0+d2+d02
local=(R00+R10+R01+R11)/4;extreme=(R00+R11)/2
assert s.expand(extreme-local-d02/4)==0
print('FAIL EXPECTED two-extreme global average: excess double-crossing coefficient delta02/4')
# Each binary grouping uses independent proper-face balancing, hence the same four terms.
first0=((R00+R10)/2+(R01+R11)/2)/2
first2=((R00+R01)/2+(R10+R11)/2)/2
assert s.expand(first0-first2)==0
print('PASS independent local forest balancing diamond')
# Source0 full mixed moment: arbitrary current/source spectators, no constant-C shortcut.
h,t,r1,r2=s.symbols('h t r1 r2')
sourcejet=-1/r1+(h-1)/r2
Ceff=-1/r1+(t-1)/r2
assert s.factor(sourcejet-Ceff)==(h-t)/r2
print('PASS Y0 spectator functional is evanescent h-t, including external current')
print('PASS total 4; one retained EXPECTED FAIL shortcut')
