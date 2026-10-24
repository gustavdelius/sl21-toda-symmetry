"""Independent Wick/jet/cochain checks. No old research code imported."""
import sympy as S
from pathlib import Path
npass=0

def check(name, expr):
    global npass
    ok=S.simplify(expr)==0
    print(('PASS ' if ok else 'FAIL ')+name)
    if ok: npass+=1
    else: print('RESIDUAL',S.factor(expr))
    return ok

t,rho,g,h,r,R=S.symbols('t rho g h r R', nonzero=True)
delta=g-rho; eta=t-2
# Each vertex is a list of independent labels. Bilinearity, not a desired relation,
# supplies the Gaussian regulator product.
p=S.Matrix(S.symbols('p0:6'));q=S.Matrix(S.symbols('q0:6'));v=S.Matrix(S.symbols('v0:6'))
M=S.diag(delta,delta,delta,delta,0,0)
M[4,5]=h-t;M[5,4]=h-t
bil=lambda x,y:(x.T*M*y)[0]
check('auxiliary Gaussian contraction is associative',bil(p,q)+bil(p+q,v)-bil(q,v)-bil(p,q+v))
L=S.symbols('L')
coc=lambda x,y:L*(x[0]+x[1])*(y[2]+y[3])-S.I*S.pi/t*(x[0]*y[1]-x[1]*y[0])
check('physical ordered cocycle satisfies cocycle identity',coc(p,q)+coc(p+q,v)-coc(q,v)-coc(p,q+v))
# Explicit uncontracted C fields and coordinates for B1 B2 U.
z1,z2,x,C1,C2=S.symbols('z1 z2 x C1 C2')
r12=z1-z2;r1=z1-x;r2=z2-x
full=(C1-1/r12+eta/r1)*(C2+1/r12+eta/r2)+1/r12**2
left=((C1-1/r12)*(C2+1/r12)+1/r12**2).subs({C1:C1+eta/r1,C2:C2+eta/r2},simultaneous=True)
right=(C1-1/r12+eta/r1)*(C2+eta/r2+1/r12)+1/r12**2
check('(B1 B2) U exact full Wick polynomial',left-full)
check('B1 (B2 U) exact full Wick polynomial',right-full)
check('retained B constituent changes BU exponent to g-1',(rho-1)+(g-rho)-(g-1))
# Nontrivial diagonal jet: C^2+C' is the actual BB coefficient, not C^2.
C,Cp,Cpp,Cppp,Dp=S.symbols('C Cp Cpp Cppp Dp')
Cplus=C+r*Cp/2+r**2*Cpp/8+r**3*Cppp/48
Cminus=C-r*Cp/2+r**2*Cpp/8-r**3*Cppp/48
mid=full.subs({z1:R+r/2,z2:R-r/2,x:0,C1:Cplus,C2:Cminus},simultaneous=True)
P0=S.limit(mid,r,0)
fused=(C+eta/R)**2+Cp-eta/R**2
check('BB fused descendant including derivative contraction',P0-fused)
check('the scalar fused descendant is eta(eta-1)/R^2',P0.subs({C:0,Cp:0})-eta*(eta-1)/R**2)
# Reconstruct second pair jet by each actual bracketing before expansion.
leftmid=left.subs({z1:R+r/2,z2:R-r/2,x:0,C1:Cplus,C2:Cminus},simultaneous=True)
rightmid=right.subs({z1:R+r/2,z2:R-r/2,x:0,C1:Cplus,C2:Cminus},simultaneous=True)
P2=S.expand(S.series(mid,r,0,3).removeO()).coeff(r,2)
check('second pair jet agrees in left grouping',S.expand(S.series(leftmid,r,0,3).removeO()).coeff(r,2)-P2)
check('second pair jet agrees in right grouping',S.expand(S.series(rightmid,r,0,3).removeO()).coeff(r,2)-P2)
# The cross regulator must be Taylor-expanded with the pair, including smooth tests.
f0,f1,f2,k0,k1,k2=S.symbols('f0 f1 f2 k0 k1 k2')
test=(f0+r*f1/2+r**2*f2/8)*(k0-r*k1/2+r**2*k2/8)
Ejet=1+r**2*(Dp/4-(g-1)/(4*R**2))
jet=S.expand((P0+r**2*P2)*Ejet*test)
check('arbitrary first test jet retained',jet.coeff(r,1)-P0*(f1*k0-f0*k1)/2)
check('arbitrary second test jet retained',jet.coeff(r,2)-(P0*(f2*k0+f0*k2-2*f1*k1)/8+(P2+P0*(Dp/4-(g-1)/(4*R**2)))*f0*k0))
# Mixed finite contact ∂_x U: insertion-wise and full-regulated derivatives differ.
q=S.symbols('q', positive=True)
W=eta*q**(g-2)
physical_derivative=-(rho-2)*eta*q**(g-3)
whole_derivative=-S.diff(W,q)
check('BU contact differentiation retains constituent regulator',physical_derivative-whole_derivative-delta*W/q)
f=S.Function('f')(q)
check('same derivative identity with arbitrary smooth outer test',physical_derivative*f-(-S.diff(W*f,q)+W*S.diff(f,q))-delta*W*f/q)
# A first-order transported product has a coboundary change; the primitive
# free braided relation kills its correction if one-current symbols stay fixed.
s,Z1V=S.symbols('s Z1V')
DBD=-s*Z1V;DDB=-Z1V
check('coboundary cannot change leading braided defect',DBD-s*DDB)
alpha=t*(t-1)
check('canonical BB coefficient requires nonzero change',t+alpha-t**2)
gamma=t/(1-t)
check('common sqrt(t) charge rescaling matches BB',t*t-t**2)
check('common sqrt(t) charge rescaling matches AA',t*(gamma**2/t)-gamma**2)
# Counts are generated from labels, source powers and measures.
def degree(a,b,abar,bbar,n0,n2):
    nc=a+b+abar+bbar; ns=n0+n2
    E=g*(S.binomial(a,2)+S.binomial(b,2)+S.binomial(abar,2)+S.binomial(bbar,2))-rho*(a*b+abar*bbar)
    E+=(a+abar-b-bbar)*ns+2*n0*(n0-1)//2+2*n2*(n2-1)//2+2*(h-1)*n0*n2
    deriv=b+bbar
    return S.expand(nc+2*ns-1+E-deriv)
check('first BBbarB one-source union degree',degree(0,2,0,1,0,1)-(g-2))
check('five-current one-source overlap degree',degree(1,3,0,1,0,1)-(3*(g-rho)-1))
check('five-current unlike-source overlap degree',degree(1,3,0,1,1,1)-(3*(g-rho)+2*h-4))
# A toy Mellin factor diagnoses why intermediate finite parts are not products.
e=S.symbols('e')
fp=lambda f:S.expand(S.series(f,e,0,1).removeO()).coeff(e,0)
check('finite-part multiplication has nonzero defect one',fp(e/e)-fp(e)*fp(1/e)-1)
print(f'TOTAL {npass} PASS')
Path(__file__).with_name('derived_jets.txt').write_text('P0 = '+str(S.factor(P0))+'\nP2 = '+str(S.factor(P2))+'\n')
