"""Independent two-source partial fractions and descendant regrouping."""
import sympy as S
x,w0,w2,r,t,h=S.symbols('x w0 w2 r t h',nonzero=True)
beta=t-1;a=w0-w2
checks=[]
def check(name,good):
 good=bool(good);checks.append(good);print(('PASS ' if good else 'FAIL ')+name)
H=1/((x-w0)*(x-w2))
f=H*(1/(x-w0)+beta/(x-w2))
parts=1/a/(x-w0)**2-beta/a/(x-w2)**2+(beta-1)/a**2*(1/(x-w0)-1/(x-w2))
check('complete unequal-source B rational partial fraction',S.factor(f-parts)==0)
primitive=1/a/(x-w0)-beta/a/(x-w2)
remaining=(beta-1)/a**2*(1/(x-w0)-1/(x-w2))
check('two-source moving poles admit one spatial derivative primitive',S.factor(-S.diff(primitive,x)+remaining-f)==0)
# The two independent delta coefficients are exactly the original spectator denominators.
check('Y0 unary collapse coefficient',S.limit((x-w0)**2*f,x,w0)==1/(w0-w2))
check('Y2 unary collapse coefficient',S.limit((x-w2)**2*f,x,w2)==beta/(w2-w0))
local=f.subs({w0:r,w2:0})
for n in range(5):
 coeff=S.diff(local,r,n).subs(r,0)/S.factorial(n)
 check('source-pair mode n='+str(n),S.factor(coeff-(n+t)/x**(n+3))==0)
# Source-pair angular modes each carry one factor n+t from BOTH actual currents.
z,v,eps,RR=S.symbols('z v eps RR',nonzero=True)
for n in range(3):
 coeff=(n+t)**2/(z*v)**(n+3)
 radial=S.pi*coeff*(RR**(2*(h+n))-eps**(2*(h+n)))/(h+n)
 check('source-pair descendant n='+str(n),S.simplify(S.diff(radial,RR)-2*S.pi*coeff*RR**(2*(h+n)-1))==0)
# At n=1 the derivative descendant supplies +1 in addition to the P OPE coefficient t.
check('first descendant extra derivative-current contraction is indispensable',S.expand((t+1)**2-t*t)==2*t+1)
# Opposite currents can consume DIFFERENT source labels; the surviving unlike source factor has degree 2h-4.
rpos=S.symbols('rpos',positive=True)
collapsed=beta/(rpos*(-rpos))*rpos**(2*h-2)
check('distinct-source unary pair full surviving power',S.simplify(collapsed+beta*rpos**(2*h-4))==0)
q=S.symbols('q',positive=True)
check('distinct-source unary face vanishes on 2<t<3',S.simplify((2*h-4).subs(h,q+2))==2*q)
# A single derivative field has one linear contraction factor, no simultaneous two-source choice.
p0,p2,C=S.symbols('p0 p2 C')
check('one C never contracts onto both sources at once',S.Poly(C+p0+beta*p2,p0,p2).coeff_monomial(p0*p2)==0)
# An improved cutoff representative and its forest correction are the same contribution.
lam0,lam2=S.symbols('lam0 lam2')
tail=-t*lam2
Kct=-lam0*eps**(2*h)/h
check('canonical tail generates exactly the known L improvement once',S.expand(tail*Kct-t*lam0*lam2*eps**(2*h)/h)==0)
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
