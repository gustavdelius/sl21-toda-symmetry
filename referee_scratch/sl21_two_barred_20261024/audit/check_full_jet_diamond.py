import sympy as s
x,y,u1,v1,u2,v2=s.symbols('x y u1 v1 u2 v2')
q1,q2=s.symbols('q1 q2')
psi=s.Function('psi')(x,y)
C1=s.Function('C1')(x,y); C2=s.Function('C2')(x,y)
D1=s.Function('D1')(x,y); D2=s.Function('D2')(x,y)
H=s.Function('H')(x,y); Hb=s.Function('Hb')(x,y)
G=psi*(1+C1*u1+D1*v1+C2*u2+D2*v2+(C1*C2+H)*u1*u2+C1*D2*u1*v2+D1*C2*v1*u2+(D1*D2+Hb)*v1*v2)
zero={u1:0,v1:0,u2:0,v2:0}
def M1(f):return q1*q1*s.diff(f,x)-q1*(s.diff(f,u1)+s.diff(f,v1))
def M2(f):return q2*q2*s.diff(f,y)-q2*(s.diff(f,u2)+s.diff(f,v2))
def ck(name,e):
 z=s.expand(e)
 assert z==0,(name,z)
 print('PASS '+name,flush=True)
a=s.expand(M1(M2(G)).subs(zero)); b=s.expand(M2(M1(G)).subs(zero))
ck('full Gaussian disjoint diamond with arbitrary source-dependent coefficients',a-b)
expected=q1**2*q2**2*s.diff(psi,x,y)-q1**2*q2*s.diff((C2+D2)*psi,x)-q1*q2**2*s.diff((C1+D1)*psi,y)+q1*q2*((C1+D1)*(C2+D2)+H+Hb)*psi
ck('complete double contact formula including both Wick crossmatchings',a-expected)
ck('source derivative of second coefficient retained',a.coeff(s.diff(C2,x))+q1**2*q2*psi)
ck('source derivative of first coefficient retained',a.coeff(s.diff(C1,y))+q1*q2**2*psi)
ck('holomorphic crossmatching coefficient',a.coeff(H)-q1*q2*psi)
ck('antiholomorphic crossmatching coefficient',a.coeff(Hb)-q1*q2*psi)
# Arbitrary firstjets of each source include a mixed secondjet of the joint test.
ck('mixed second source jet coefficient',a.coeff(s.diff(psi,x,y))-q1**2*q2**2)
# Cocycle associativity: exponent is bilinear in full constituent labels.
ns=s.symbols('n0:4'); ms=s.symbols('m0:4'); ks=s.symbols('k0:4')
add=lambda a,b:tuple(x+y for x,y in zip(a,b))
E=lambda a,b:(a[0]+a[1])*(b[2]+b[3])
F=lambda a,b:a[0]*b[1]-a[1]*b[0]
for label,Q in [('sign',E),('phase',F)]:ck('full constituent cocycle '+label+' associative',Q(ns,ms)+Q(add(ns,ms),ks)-Q(ms,ks)-Q(ns,add(ms,ks)))
# Independent product average versus global two-extreme average.
R,d1,d2,d12=s.symbols('R d1 d2 d12')
f00=R;f10=R+d1;f01=R+d2;f11=R+d1+d2+d12
ck('extreme average misses quarter of disjoint double crossing',(f00+f11)/2-(f00+f10+f01+f11)/4-d12/4)
print('SUMMARY 10 PASS, 0 FAIL')
print('FULL MOMENT:',expected)
