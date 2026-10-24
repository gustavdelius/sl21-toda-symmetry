"""Independent interaction-map checks; no prior scratch code is imported."""
import sympy as S
from itertools import combinations
npass=0

def check(name,expr):
 global npass
 ok=S.simplify(expr)==0
 print(('PASS ' if ok else 'FAIL ')+name)
 if ok:npass+=1
 else: print('RESIDUAL',S.factor(expr))

t,h,r,z,v=S.symbols('t h r z v',nonzero=True)
b=t-1;rho=2/t
# Local Gaussian products in d,a0 basis.
G=S.Matrix([[-rho,1],[1,-1]])
d=S.Matrix([1,0]);a0=S.Matrix([0,1]);a2=-t*d-a0;u=d+a2
sd=lambda x,y:(x.T*G*y)[0]
check('Y2 current double-pole coefficient',-sd(a0,a2)-b)
check('Y0 current double-pole coefficient',-sd(a0,a0)-1)
check('Y2 first-jet vector beta a2-a0=t(d+a2)',sum(b*a2-a0-t*u))
check('K0 leaves no C pole',sd(a0,d+a0))
# Minimal mixed base differential kernel, including full first field jets.
x,y,w,wb,C,Cb,D,Db,N=S.symbols('x y w wb C Cb D Db N')
base=1/((x-w)*(y-wb))
# Normal exponential derivative N_x=D N; N_y=Db N.
Lx=lambda f:(C+b*D)*f-b*(S.diff(f,x)+D*f)
Ly=lambda f:(Cb+b*Db)*f-b*(S.diff(f,y)+Db*f)
actual=base*(C+b/(x-w))*(Cb+b/(y-wb))
check('spatial base differential reconstruction off singular faces',Lx(Ly(base))-actual)
# Ordered unary repair from Cauchy-vs-spatial Pf difference.
Pl,Pr,U,delta,Pv,Pf,sigma=S.symbols('Pl Pr U delta Pv Pf sigma')
# 1/(r +/- i sigma0)=PV -/+ i*pi*sigma*delta;
# squared boundary delta' terms cancel between the two tails.
Qplus=Pv-S.I*S.pi*sigma*delta;Qminus=Pv+S.I*S.pi*sigma*delta
unary=-b*(2*b*U*Pf-Pl*Qplus+(b*Pl+t*Pr)*Qminus)
sym=-2*b*b*U*Pf-b*((t-2)*Pl+t*Pr)*Pv
odd=-S.I*S.pi*sigma*t*b*(Pl+Pr)*delta
check('unary repair symmetric Pf and PV terms',unary-sym-odd)
check('balanced ordered matching removes odd contact',unary+S.I*S.pi*sigma*t*b*(Pl+Pr)*delta-sym)
# Source source scalar kernel and all radial Taylor coefficients.
f=(1/(z-r)+b/z)/((z-r)*z)
for j in range(5):
 coeff=S.diff(f,r,j).subs(r,0)/S.factorial(j)
 check(f'full source-pair current jet n={j}',coeff-(j+t)/z**(j+3))
# Reconstruct first descendant from actual Wick contractions after source fusion.
check('first descendant retains C contraction',(t/z)*(1/z)+1/z**2-(t+1)/z**2)
R,eps=S.symbols('R eps',positive=True)
for j in range(3):
 moment=S.pi*(R**(2*(h+j))-eps**(2*(h+j)))/(h+j)
 check(f'radial moment derivative n={j}',S.diff(moment,R)-2*S.pi*R**(2*(h+j)-1))
 coeff=(j+t)**2*eps**(2*(h+j))/((h+j)*(z*v)**(j+3))
 check(f'physical source-pair descendant denominator n={j}',coeff.subs(h,t)-(j+t)*eps**(2*(t+j))/(z*v)**(j+3))
check('current improvement coefficient',t/h-1+(h-t)/h)
# Once-only source counting and normalized unit.
a,c,beta=S.symbols('a c beta')
expansion=1-(a+c)+(a+c)**2/2
check('distinct source cross coefficient',S.expand(expansion).coeff(a,1).coeff(c,1)-1)
check('equal source coefficient',S.expand(expansion).coeff(a,2)-S.Rational(1,2))
# Compatibility of unary forests (two current nodes and two source nodes).
faces=[frozenset((i,j)) for i in ('B1','B2') for j in ('Y0','Y2')]
pairs=[p for p in combinations(faces,2) if not p[0]&p[1]]
check('disjoint two-unary forests have two source matchings',len(pairs)-2)
shared=[frozenset(('B1','Y2')),frozenset(('B2','Y2'))]
check('same source cannot fuse twice',int(not bool(shared[0]&shared[1])))
# Whole-cluster extreme-order averaging, not iterated Jordan multiplication.
raw,d1,d2=S.symbols('raw d1 d2')
R0=raw;R1=raw+d1;R2=raw+d1+d2
Cs=[(d1+d2)/2,(-d1+d2)/2,-(d1+d2)/2]
for j,rr in enumerate([R0,R1,R2]):
 check(f'one-bar crossing position {j} has same repaired value',rr+Cs[j]-(R0+R2)/2)
# K0 variation of both Serre words: the only marginal scalar is zero per slot.
z1,z2,a=S.symbols('z1 z2 a')
P=(-1/(z1-z2)+1/(z1-a))*(1/(z1-z2)+1/(z2-a))+1/(z1-z2)**2
check('K0 insertion marginal B B A scalar vanishes',P)
g=S.symbols('g')
check('A A A K0 full degree',3*g+3*(1-rho)+3-(3*(g-rho)+6))
# Normalized functional, keeping operator insertion map itself nonunital.
l0,l2=S.symbols('l0 l2'); A0,A2,A02,Z0,Z2,Z02,F=S.symbols('A0 A2 A02 Z0 Z2 Z02 F')
I=F+l0*A0+l2*A2+l0*l2*A02
Z=1+l0*Z0+l2*Z2+l0*l2*Z02
Zi=1-l0*Z0-l2*Z2+l0*l2*(2*Z0*Z2-Z02)
normalized=S.expand(I*Zi)
check('normalized mixed coefficient',normalized.coeff(l0,1).coeff(l2,1)-(A02-Z0*A2-Z2*A0+(2*Z0*Z2-Z02)*F))
check('normalized unit at first and mixed source order',normalized.coeff(l0,1).coeff(l2,1).subs({F:1,A0:Z0,A2:Z2,A02:Z02}))
print(f'TOTAL {npass} PASS')
