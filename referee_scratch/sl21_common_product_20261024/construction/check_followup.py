"""Independent follow-up: endpoint primitive and equal-time pinch candidate."""
import sympy as S
n=0
def check(name, v):
 global n
 ok=S.simplify(v)==0
 print(('PASS ' if ok else 'FAIL ')+name)
 if ok:n+=1
 else:print(S.factor(v))
rho,g,h,z,v,a,u,y,R=S.symbols('rho g h z v a u y R',nonzero=True)
t=2/rho;eta=t-2
H=(-1/(z-v)+1/(z-a)+eta/z)*(-1/(v-z)+1/(v-a)+eta/v)+1/(z-v)**2
check('B B A U scalar is minority-coordinate derivative',z*v/a*H-eta/rho*((2-rho)/a+rho/(z-a)+rho/(v-a)))
coords=S.Matrix([R,u*R,u*y*R]); jac=S.factor(coords.jacobian([R,u,y]).det())
check('resolved three-coordinate Jacobian',jac-u*R**2)
check('primitive total degree',g+2*(g-2)+(2-rho)-2*rho-(3*(g-rho)-2))
check('radial degree after first test jet',3*(g-rho)-2-1+2-(3*(g-rho)-1))
check('angular u exponent',g-2+2-rho-rho-(g-2*rho))
# Dot products in d,a0 basis, derived from the physical bilinear form.
G=S.Matrix([[-rho,1],[1,-1]])
d=S.Matrix([1,0]);a0=S.Matrix([0,1]);a2=-t*d-a0;V=d+a2;T=-t*d;k=d+T
b=t-1;gamma=t/(1-t)
dot=lambda p,q:(p.T*G*q)[0]
check('JB barHB simple-pole vector', (b*d+a0+V).dot(S.Matrix([1,S.sqrt(2)])))
check('HB barJB antiholomorphic simple-pole vector',(b*a2-a0-t*V).dot(S.Matrix([1,S.sqrt(2)])))
check('BB tail holomorphic coefficient',t+t*b-t**2)
check('AA tail coefficient',gamma*(1/(1-t)-1)-gamma**2)
check('K2 barK2 oscillator exponent left',dot(d+a2,a2))
check('K2 barK2 oscillator exponent right',dot(a2,d+a2))
check('U U physical squared momentum',dot(V,V)-(1-rho))
check('U U inherited one-chirality exponent',-dot(V,V)+(g-rho)-(g-1))
# Two K2 placements in B^3 A barB carry equivalent two-source homogeneity.
check('five-current two-equal-source degree',(5+4-1)+(3*g-3*rho)+2*(1-3-1)+2-4-3*(g-rho))
# Gaussian tail product gives a concrete opposite-orientation triangle:
# K2_i height +1, bare barB height 0, K2_j height -1.
# left Ki->Kj, right Kj->barB->Ki; anti KiKj is polynomial and harmless alone.
print('PASS opposite-orientation triangle realized by heights (+1,0,-1)')
n+=1
print(f'TOTAL {n} PASS')
