"""Fresh local Ward normalization and even-root collision calculations."""
import sympy as S
import numpy as np
from scipy.integrate import quad
from scipy.special import iv

passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    passed+=int(bool(ok));failed+=int(not bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)
def eq(name,a,b=0):check(name,S.expand(S.together(a-b).as_numer_denom()[0])==0)

t,h,rho=S.symbols('t h rho',nonzero=True)
C,D,C2,D2=S.symbols('C D C2 D2')
# The chiral source point splittings are expanded, not copied as N.
A2=-t*D-C;A22=-t*D2-C2
Q0=(A22+A2*A2)/2;Q2=(C2+C*C)/2
N=S.expand(Q2-Q0)
K=(1-t)*D;K2=(1-t)*D2
pole3=t;pole2=C+t*D;pole1=C2+C*D+t*(D2+D*D)/2
response=-pole1+(C2+t*D2+pole2*K)-t*(K2+K*K)/2
RP=-t*C*D+t*t*D2/2-t**3*D*D/2
eq('triple/double/simple potential poles give full distributional response',response,RP)
eq('double-Yukawa contact plus potential equals derivative of C E',N-RP/t,C2+(1-t)*C*D)
# F_A has momentum -k, whereas E above has momentum +k: dF_A/F_A=-K.
eq('even-root A Ward coefficient',(-K)-(-D),(-t/(1-t))*(-K))

# Independent analytic-regulator subtractions, with ordinary area and e^-S.
R,eps=S.symbols('R eps',positive=True)
radial=S.pi*(R**(2*h)-eps**(2*h))/h
eq('source-pair action subtraction coefficient',radial/S.pi**2+eps**(2*h)/(S.pi*h),R**(2*h)/(S.pi*h))
eq('tail insertion subtraction coefficient',t*radial/S.pi+t*eps**(2*h)/h,t*R**(2*h)/h)
eq('physical current subtraction is recovered only after h=t',(t/h).subs(h,t),1)

# Circle derivative of the exclusion selects the quadratic holomorphic jet.
tt=np.sqrt(5.);cv=.23+.11j;dv=.17-.08j;c2v=.07+.03j;d2v=-.04+.05j
nval=complex(N.subs({t:tt,C:cv,D:dv,C2:c2v,D2:d2v}))
theta=2*np.pi*(np.arange(8192)+.5)/8192;unit=np.exp(1j*theta)
errors=[]
for ep in (.04,.02,.01):
    rr=ep*unit
    G2=np.exp(cv*rr+c2v*rr*rr/2+(.09+.02j)*np.conj(rr))
    G0=np.exp((-tt*dv-cv)*rr+(-tt*d2v-c2v)*rr*rr/2+(.05-.04j)*np.conj(rr))
    value=np.mean(np.exp(-2j*theta)*(G2-G0))/ep**2
    errors.append(abs(value-nval))
    print(f'DATA Ward circle epsilon={ep} normalized_error={errors[-1]:.8g}',flush=True)
check('direct regulated Ward circle converges to N',errors[-1]<errors[0]/12)
check('Ward circle remainder has two additional powers',abs(errors[-1]/errors[-2]-.25)<.01)

# Negative-h tests use nonconstant source backgrounds and an independently
# integrated radial Bessel function. The matched counterterms fix the finite
# parts; the overall reference length is one.
hh=-np.sqrt(2)/5;tt=np.sqrt(5.);ar=.37;br=.23;outer=.8
target_action=sum(np.pi*outer**(2*(hh+n))/((hh+n)*float(S.factorial(n))**2) for n in range(30))
target_tail=sum(np.pi*ar**(n+1)*br**n*outer**(2*(hh+n))/((hh+n)*float(S.factorial(n+1))*float(S.factorial(n))) for n in range(30))
action_errors=[];tail_errors=[]
for ep in (.08,.04,.02,.01):
    raw=2*np.pi*quad(lambda x:x**(2*hh-1)*iv(0,2*x),ep,outer,epsabs=1e-10)[0]
    action=raw+np.pi*ep**(2*hh)/hh
    # e^(a*r+b*rbar) has the selected m=n+1 moment sqrt(a/b) I1.
    raw_tail=2*np.pi*np.sqrt(ar/br)*quad(lambda x:x**(2*hh-2)*iv(1,2*np.sqrt(ar*br)*x),ep,outer,epsabs=1e-10)[0]
    tail=tt/np.pi*raw_tail+tt*ar*ep**(2*hh)/hh
    action_errors.append(abs(action-target_action));tail_errors.append(abs(tail-tt/np.pi*target_tail))
    print(f'DATA subtraction epsilon={ep} action_error={action_errors[-1]:.8g} tail_error={tail_errors[-1]:.8g}',flush=True)
for name,errors in [('action',action_errors),('tail',tail_errors)]:
    check(name+' matched subtraction converges',errors[-1]<errors[0]/10)
    check(name+' next descendant has power 2h+2',abs(errors[-1]/errors[-2]-2**(-(2*hh+2)))<.01)

# The new potentially marginal channel is AAAB/V1, not the unlike double
# Yukawa primitive at generic t. Reconstruct it at general minority B point z.
z,w=S.symbols('z w');apos=[S.Rational(-5,7),S.Rational(1,9),S.Rational(8,7)]
rr=2/t;logR=-rr*sum(1/(z-a) for a in apos)
G=(z-w)**2/S.prod((a-w)**2 for a in apos)
source=G*(sum(1/(z-a) for a in apos)-t/(z-w))
eq('even-root scalar kernel is a B-coordinate derivative',source,-(S.diff(G,z)+logR*G)/rr)
for i,ai in enumerate(apos):
    Di=1/S.prod((ai-aj)**2 for j,aj in enumerate(apos) if j!=i)
    Si=(z-ai)**2*Di
    f=sum(1/(z-aj) for aj in apos)-t/(z-ai)
    eq(f'even-root tail derivative at A slot {i}',Si*f,-(S.diff(Si,z)+logR*Si)/rr)
gamma=t/(1-t)
eq('combined even-root scalar primitive normalization',1+gamma,1/(1-t))

# Isolated AA pair/triple coefficients. F has no derivative polynomial;
# compute its plane finite part from actual residues at distinct real points.
for count,pos in [(2,[S.Rational(-2,5),S.Rational(3,5)]),
                  (3,[S.Rational(-4,5),S.Rational(1,7),S.Rational(6,5)])]:
    area=0;tail=0
    for i,zi in enumerate(pos):
        Q=1/S.prod((w-zj)**2 for j,zj in enumerate(pos) if j!=i)
        area-=zi*S.diff(Q,w).subs(w,zi)
        tail+=Q.subs(w,zi)
    eq(f'AA-type source area n={count}',area,tail)
    eq(f'AA-type action and all tails n={count}',-area-gamma*tail,-tail/(1-t))
r,Lsep=S.symbols('r Lsep',nonzero=True)
pos=[r/2,-r/2,-Lsep]
summed=sum(1/S.prod((zi-zj)**2 for j,zj in enumerate(pos) if j!=i) for i,zi in enumerate(pos))
local=S.cancel(-r*r*summed/(1-t))
eq('new AA pair embedded in AAA has the same coefficient',local.subs(r,0),-2/((1-t)*Lsep**2))
eq('new AA pair embedding has no odd midpoint pole',S.diff(local,r).subs(r,0))
# The two-majority-pole primitive has a nonzero infinity moment.
pa=[S.Rational(-2,3),S.Rational(4,5)]
area=0;Ssum=0
for i,ai in enumerate(pa):
    Q=(z-w)**2/S.prod((w-aj)**2 for j,aj in enumerate(pa) if j!=i)
    area-=ai*S.diff(Q,w).subs(w,ai)
    Ssum+=Q.subs(w,ai)
eq('proper AAB scalar area retains the infinity constant',area,Ssum-1)
eq('proper AAB combined source/tail primitive',area+gamma*Ssum,(Ssum+t-1)/(1-t))
print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
