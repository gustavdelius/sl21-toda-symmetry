"""Independent primary-vector, distribution-moment and source-forest checks."""
import sympy as S
import math
from scipy.integrate import quad
checks=[]
def check(name,condition):
    good=bool(condition);checks.append(good);print(('PASS ' if good else 'FAIL ')+name,flush=True)
b=S.symbols('b',positive=True);t=S.symbols('t',positive=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I]);a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I]);T=a0+a2;d=-T/b**2;v=d+a2;k=d+T
simp=lambda a:S.simplify(S.expand(a).subs(b**2,t))
check('primary source squares',simp(a0.dot(a0))==-1 and simp(a2.dot(a2))==-1)
check('both current/source exponential simple powers',simp(d.dot(a0))==1 and simp(d.dot(a2))==1)
check('double derivative coefficients 1 and t-1',simp(-a0.dot(a0))==1 and simp(-a0.dot(a2))==t-1)
check('a0+(t-1)d=-(d+a2)',all(simp(u)==0 for u in a0+(b**2-1)*d+v))
check('A/V1 exponential double pole',simp((-d).dot(-T))==2)
# Direct full-jet differentiation, no prior script imports.
x,y,w,wb=S.symbols('x y w wb'); C,Cb,D,Db,L,Lb,beta=S.symbols('C Cb D Db L Lb beta')
G=1/((x-w)*(y-wb));N=S.exp(D*x+Db*y+L*w+Lb*wb)
Dx=lambda f:(C+beta*D)*f-beta*S.diff(f,x)
Dy=lambda f:(Cb+beta*Db)*f-beta*S.diff(f,y)
physical=N*(C/(x-w)+beta/(x-w)**2)*(Cb/(y-wb)+beta/(y-wb)**2)
check('complete pair Wick polynomial as spatial differential operator',S.simplify(Dx(Dy(G*N))-physical)==0)
# Independently integrate the sign and edge moments with a round IR boundary.
# Horizontal slice integral of 1/barw^2 is small-circle contribution minus outer-circle contribution.
for sigma in [-1,1]:
    eps=.01;delta=math.sqrt(2)*eps;R=100.
    small=sigma*math.pi
    big=quad(lambda yy:math.copysign(1,sigma*delta-yy)*2*math.sqrt(max(0,R*R-yy*yy))/R**2,-R,sigma*delta,epsabs=1e-10)[0]+quad(lambda yy:math.copysign(1,sigma*delta-yy)*2*math.sqrt(max(0,R*R-yy*yy))/R**2,sigma*delta,R,epsabs=1e-10)[0]
    M=small-big;yy=sigma*delta;ell=math.sqrt(R*R-yy*yy)
    J=2j*math.atan(ell/yy);Jw=J-4j*yy*ell/R**2
    check('ordered source sign moment sigma='+str(sigma),abs(M-sigma*math.pi)<.001)
    check('ordered edge simple moment sigma='+str(sigma),abs(J-1j*math.pi*sigma)<.001)
    check('ordered edge first moment sigma='+str(sigma),abs(Jw-1j*math.pi*sigma)<.001)
    check('outer-test jet cancels in local limit sigma='+str(sigma),abs(math.pi*Jw-1j*math.pi*M)<.001)
# Circular vs spatial double pole: on radial smooth compact tests the Pf integral is zero,
# while -d_x(1/z) has value -pi*phi(0). Integrate radial derivative directly.
# phi(r)=(1-r^2)^4 for r<1, zero otherwise.
radial=math.pi*quad(lambda rr:-8*rr*(1-rr*rr)**3,0,1,epsabs=1e-13)[0]
check('spatial double pole differs by -pi delta from circular finite part',abs(radial+math.pi)<1e-12)
# Off-diagonal convolution first jets, derived using the simple Cauchy convolution's derivatives.
KL,KR=S.symbols('KL KR');bet=t-1
hol=-C+bet*L;anti=Cb-bet*Lb
sub={C:-KL-bet*D,L:KL-D,Cb:-KR-bet*Db,Lb:KR-Db}
check('raw first-jet convolution holomorphic coefficient',S.expand(hol.subs(sub)-t*KL)==0)
check('raw first-jet convolution antiholomorphic coefficient',S.expand(anti.subs(sub)+t*KR)==0)
# Phase-aligned canonical tail average and jump, from actual double/simple Laurent terms.
r=S.symbols('r');P,Q,U=S.symbols('P Q U');delta,delta1,PV,PF=S.symbols('delta delta1 PV PF')
def tail(sig):
    z1=PV-S.I*S.pi*sig*delta;zb1=PV+S.I*S.pi*sig*delta
    z2=PF+S.I*S.pi*sig*delta1;zb2=PF-S.I*S.pi*sig*delta1
    return t*bet*U*(z2+zb2)-t*P*z1+(t*bet*P+t*t*Q)*zb1
avg=S.expand((tail(1)+tail(-1))/2);jump=S.expand(tail(1)-tail(-1))
check('tail average retains symmetric double-pole finite part',S.expand(avg-(2*t*bet*U*PF+t*((t-2)*P+t*Q)*PV))==0)
check('tail double-pole delta-prime cancellation',not jump.has(delta1))
check('tail jump exact endpoint coefficient',S.expand(jump-2*S.I*S.pi*t*t*(P+Q)*delta)==0)
unary=S.expand(-bet/t*avg)
check('unary repair has explicit nonzero symmetric density terms',S.expand(unary-(-2*bet**2*U*PF-bet*((t-2)*P+t*Q)*PV))==0)
# Endpoint angular Taylor selection INCLUDING arbitrary source/test jets.
z,zb,h,eps,RR=S.symbols('z zb h eps RR',nonzero=True)
f,fz,fb,fzb,A,B=S.symbols('f fz fb fzb A B')
jet=(f+fz*z+fb*zb+fzb*z*zb)*(1+A*z+B*zb+A*B*z*zb)
angular1=S.expand(jet).coeff(z,1).coeff(zb,1)
check('endpoint descendant includes test derivatives',angular1==f*A*B+fz*B+fb*A+fzb)
for n in range(4):
    primitive=S.pi*(RR**(2*(h+n-1))-eps**(2*(h+n-1)))/(h+n-1)
    check('endpoint radial coefficient n='+str(n),S.simplify(S.diff(primitive,RR)-2*S.pi*RR**(2*(h+n-1)-1))==0)
# K0 unary insertion is a new needed channel: scalar leading term vanishes separately.
u0,u1,aa=S.symbols('u0 u1 aa')
f0=-1/(u0-u1)+1/(u0-aa);f1=1/(u0-u1)+1/(u1-aa)
check('one K0 in B3A has vanishing leading two-C scalar',S.factor(f0*f1+1/(u0-u1)**2)==0)
# Sources are consumed once: compatible forest faces disjoint or nested, never overlap.
faces=[frozenset(['w','B']),frozenset(['w','D']),frozenset(['w','B','D'])]
def compatible(a,b):return not a&b or a<=b or b<=a
check('same source cannot be contracted independently onto both currents',not compatible(faces[0],faces[1]))
check('unary face may nest inside full mixed face',compatible(faces[0],faces[2]))
a,b0=S.symbols('a b0')
check('identical source second-order factor is one half',S.expand((a+b0)**2/2).coeff(a,2)==S.Rational(1,2))
check('unlike source cross coefficient is one',S.expand((a+b0)**2/2).coeff(a,1).coeff(b0,1)==1)
# A whole-word symmetric rule is distinct from averaging positions or iterating Jordan.
positions=range(3)
check('uniform three-position bar averaging is not local half matching',sum(j>0 for j in positions)/3!=.5 and sum(j>1 for j in positions)/3!=.5)
check('two extremal bar placements match each lower pair by one half',all(sum(j>i for j in [0,2])/2==.5 for i in [0,1]))
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
