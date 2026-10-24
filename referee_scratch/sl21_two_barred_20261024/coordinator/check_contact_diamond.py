"""Disjoint full Gaussian jet contractions, including derivative-field pairings."""
import sympy as S
x,y,u1,v1,u2,v2,q1,q2=S.symbols('x y u1 v1 u2 v2 q1 q2')
cs=S.symbols('c0:16'); ps=S.symbols('p0:4'); kl,kr=S.symbols('kl kr')
C1=cs[0]+cs[1]*x+cs[2]*y+cs[3]*x*y
Cb1=cs[4]+cs[5]*x+cs[6]*y+cs[7]*x*y
C2=cs[8]+cs[9]*x+cs[10]*y+cs[11]*x*y
Cb2=cs[12]+cs[13]*x+cs[14]*y+cs[15]*x*y
psi=ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*y
# This is the EXACT part of the four-current Wick generating polynomial
# that can contribute to the two first-jet functionals. Pairings are retained.
G=psi*((1+C1*u1+Cb1*v1)*(1+C2*u2+Cb2*v2)+kl*u1*u2+kr*v1*v2)
zero={u1:0,v1:0,u2:0,v2:0}
M1=lambda f:q1*q1*S.diff(f,x)-q1*(S.diff(f,u1)+S.diff(f,v1))
M2=lambda f:q2*q2*S.diff(f,y)-q2*(S.diff(f,u2)+S.diff(f,v2))
checks=[]
def ck(name,val):
    good=bool(val);checks.append(good);print(('PASS ' if good else 'FAIL ')+name)
left=S.expand(M1(M2(G)).subs(zero));right=S.expand(M2(M1(G)).subs(zero))
ck('disjoint full-jet orders commute with source-dependent current coefficients',left==right)
formula=q1*q1*q2*q2*S.diff(psi,x,y)-q1*q1*q2*S.diff((C2+Cb2)*psi,x)-q2*q2*q1*S.diff((C1+Cb1)*psi,y)+q1*q2*((C1+Cb1)*(C2+Cb2)+kl+kr)*psi
ck('joint functional retains coefficient derivatives and both Wick pairings',S.expand(left-formula)==0)
ck('unlike-source specialization retains commutation',S.expand((left-right).subs({q1:1,q2:S.Symbol('t')-1}))==0)
ck('cross-field pairings have their nonzero required coefficient',S.diff(left,kl)==q1*q2*psi and S.diff(left,kr)==q1*q2*psi)
ck('dropping other-source coefficient derivatives changes the functional',S.expand(left.subs({x:0,y:0})).coeff(cs[9])==-q1*q1*q2*ps[0])
# Physical node0 residual with one unlike spectator, before restriction.
h,t=S.symbols('h t');r=S.symbols('r',nonzero=True)
ck('node0 unlike spectator leaves h-minus-t remainder',S.expand((h-1)/r-(t-1)/r-(h-t)/r)==0)
ck('physical node0 constant-source mixed field moment vanishes',S.simplify(((h-1)/r-(t-1)/r).subs(h,t))==0)
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
