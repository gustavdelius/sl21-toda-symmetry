"""Independent first-Yukawa boundary and descendant calculations.

Uses Gaussian Wick contractions, residue moment sums, and direct circular
boundary quadrature. Does not import any previous investigation script.
"""
import sympy as S
import numpy as np
from scipy.integrate import quad
from scipy.special import beta as betaf

passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    passed+=int(bool(ok)); failed+=int(not bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)
def eq(name,a,b=0):
    check(name,S.expand(S.together(a-b).as_numer_denom()[0])==0)

t,rho=S.symbols('t rho',nonzero=True)
eta=t-2; beta=t-1
b=S.symbols('b',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I])
a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I])
d=-(a0+a2)/b**2
dot=lambda u,v:S.simplify(u.dot(v)).subs(b**2,t)
eq('derivative self contraction',-dot(a0,a0),1)
eq('derivative/exponential contraction',dot(a0,d),1)
eq('node 2 source contraction',-dot(a0,a2),beta)
eq('node 2 fused-tail contraction',-dot(a0,d+a2),eta)

z=S.symbols('z1:4'); a,w=S.symbols('a w')
Di=[]; Ri=[]; Ti=[]
for i in range(3):
    dif=[z[i]-z[j] for j in range(3) if j!=i]
    D=1/S.prod(q*q for q in dif); l=sum(1/q for q in dif); m=sum(1/q**2 for q in dif)
    v=a-z[i]
    R=-2*D*(v+v*v*l)
    T=D*(1+4*v*l+v*v*(2*l*l+m))
    Q=(a-w)**2/S.prod((w-z[j])**2 for j in range(3) if j!=i)
    eq(f'source simple residue {i}',S.diff(Q,w).subs(w,z[i]),R)
    eq(f'source regular coefficient {i}',S.diff(Q,w,2).subs(w,z[i])/2,T)
    Di.append(D); Ri.append(R); Ti.append(T)
eq('source residue sum',sum(Ri))
eq('source first moment',sum(z[i]*Ri[i] for i in range(3)),-sum((a-z[i])**2*Di[i] for i in range(3)))

# An independent circular-boundary integration of the actual rational function.
zz=np.array([-.7+.22j,.1-.09j,.8+.04j]); aa=1.7+.31j
rad=np.array([.025,.04,.035]); arad=.03
theta=2*np.pi*np.arange(4096)/4096
unit=np.exp(1j*theta)
def G(ww):return (aa-ww)**2/np.prod((ww[:,None]-zz[None,:])**2,axis=1)
boundary=0j
for zi,rr in zip(zz,rad):
    ww=zi+rr*unit
    boundary+=-np.pi*rr*np.mean(np.conj(ww)*unit*G(ww))
ww=aa+arad*unit
A_hole=-np.pi*arad*np.mean(np.conj(ww)*unit*G(ww))
subs={a:aa,**dict(zip(z,zz))}
res=np.array([complex(R.subs(subs)) for R in Ri])
reg=np.array([complex(T.subs(subs)) for T in Ti])
expected=-np.pi*np.dot(np.conj(zz),res)-np.pi*np.dot(rad**2,reg)
check('direct circular boundary integral including holes',abs(boundary-expected)<3e-10,
      f'error={abs(boundary-expected):.3g}')
check('moving A hole contributes zero to scalar primitive',abs(A_hole)<2e-13)

# Endpoint coefficients follow by reordering A past k B insertions.
s,x=S.symbols('s x',nonzero=True)
cs=1+s+1/s
p=1-cs*x+cs*x*x-x**3
eq('endpoint polynomial factorization',p,(1-x)*(1-s*x)*(1-x/s))
eq('left endpoint cancellation',p.subs(x,s))
eq('right endpoint cancellation',p.subs(x,1/s))

# Residue moment projection. For low polynomial moments of F(w), the plane
# integral / pi equals sum_i D_i P(z_i,z_i). This is used only within its
# proved degree bound m+n < 2*N-2.
for N,pos in [(2,[S.Rational(-1,2),S.Rational(1,2)]),
              (3,[S.Rational(2,3),S.Rational(1,6),S.Rational(-5,6)])]:
    DD=[1/S.prod((pos[i]-pos[j])**2 for j in range(N) if i!=j) for i in range(N)]
    for m in range(N):
        for n in range(N-m):
            f=w**m/S.prod((w-v)**2 for v in pos)
            raw=0
            for i,zi in enumerate(pos):
                regular=w**m/S.prod((w-pos[j])**2 for j in range(N) if j!=i)
                residue=S.diff(regular,w).subs(w,zi)
                raw-=zi**(n+1)*residue/S.Integer(n+1)
            eq(f'area moment N={N}, (m,n)=({m},{n})',raw,sum(DD[i]*pos[i]**(m+n) for i in range(N)))

# Independent Wick polynomial, evaluated at finite rational configurations.
q=S.symbols('q')
C,D,C2,D2,U,U2=S.symbols('C D C2 D2 U U2')
wb=S.symbols('wb')
def wick(pos,ext):
    f=[q*C+q*q*zi*C2+sum(v/(zi-y) for y,v in ext)
       -sum(1/(zi-zj) for j,zj in enumerate(pos) if i!=j)
       for i,zi in enumerate(pos)]
    if len(pos)==1:return f[0]
    if len(pos)==2:return f[0]*f[1]+1/(pos[0]-pos[1])**2
    return f[0]*f[1]*f[2]+sum(f[k]/(pos[i]-pos[j])**2 for i,j,k in [(0,1,2),(0,2,1),(1,2,0)])
def expjets(M1,M2):return 1+q*M1+q*q*(M1*M1+M2)/2
def coefficient(f,k):return S.expand(f).coeff(q,k)

for N,pos in [(2,[S.Rational(-1,2),S.Rational(1,2)]),
              (3,[S.Rational(2,3),S.Rational(1,6),S.Rational(-5,6)])]:
    maxjet=N-1
    DD=[1/S.prod((pos[i]-pos[j])**2 for j in range(N) if i!=j) for i in range(N)]
    p2=sum(zi*zi for zi in pos)
    for node,beta_node in [(0,S.Integer(1)),(2,beta)]:
        A=C if node==0 else -t*D-C
        A2=C2 if node==0 else -t*D2-C2
        W=wick(pos,[(w,beta_node)])
        # Clear the source poles before expanding normal-field jets.
        vertex=1/S.prod(zi-w for zi in pos)
        polynomial=S.cancel(W*vertex*S.prod((w-zi)**2 for zi in pos))
        exponential=expjets(A*w+U*wb,D2*p2+A2*w*w+U2*wb*wb)
        action=[]
        for k in range(maxjet+1):
            num=coefficient(polynomial*exponential,k)
            action.append((-1)**(N+1)*sum(DD[i]*num.subs({w:zi,wb:zi}) for i,zi in enumerate(pos)))
        tails=[]
        for k in range(maxjet+1):
            tail=0
            if node==2:
                for i,zi in enumerate(pos):
                    rest=[zj for j,zj in enumerate(pos) if j!=i]
                    WW=wick(rest,[(zi,eta)])
                    EE=expjets((A+U)*zi,D2*p2+(A2+U2)*zi*zi)
                    tail+=(-1)**N*t*coefficient(WW*EE,k)/S.prod(zj-zi for zj in rest)
            tails.append(tail)
        for k in range(maxjet+1):
            actual=action[k]+tails[k]
            if node==0:expect=0
            elif N==2:expect=2*eta if k==0 else 0
            else:
                M0=sum(DD); M1=sum(zi*DD[i] for i,zi in enumerate(pos)); M2=sum(zi*zi*DD[i] for i,zi in enumerate(pos))
                F1=-t*D-C+U; F2=-t*D2-C2+U2
                expect=[-eta*(eta-1)*M0,
                        -eta*((eta-1)*F1-3*C)*M1,
                        -eta*(((eta-1)*D2/2+C2)*p2*M0+((eta-1)*(F2+F1*F1)/2-3*C*F1)*M2)][k]
            eq(f'full action+tail jet N={N}, node={node}, level={k}',actual,expect)

# Complete angular moments for the triple subtraction, including endpoints.
u=S.symbols('u')
pos=[(1+u)/3,(1-2*u)/3,-(2-u)/3]
DD=[1/u**2,1/(u*u*(1-u)**2),1/(1-u)**2]
M=[sum(pos[i]**k*DD[i] for i in range(3)) for k in range(3)]
ell=(2-u)/3; right=(1+u)/3; p2=sum(v*v for v in pos)
def angular_ratio(f):
    poly=S.Poly(S.cancel(f*u*u*(1-u)**2),u)
    return sum(co*S.rf(rho-1,k[0])/S.rf(2*rho-2,k[0]) for k,co in poly.terms())
angular={
 'I0':(M[0],(3*rho-1)/(2*rho-1)),
 'J0':(M[0]*right**2,(9*rho**2+2*rho-1)/(6*(2*rho-1)*(2*rho+1))),
 'I1r':(M[1]*right,-(3*rho+1)/(12*(2*rho-1)*(2*rho+1))),
 'IS2':(p2*M[0],rho*(3*rho+1)/((2*rho-1)*(2*rho+1))),
 'IM2':(M[2],rho*(3*rho+1)/(6*(2*rho-1)*(2*rho+1))),
}
for name,(fun,answer) in angular.items():
    eq('triple angular moment '+name,angular_ratio(fun),answer)
    rr=np.sqrt(2)
    ff=S.lambdify(u,S.cancel(fun*u*u*(1-u)**2),'numpy')
    val=quad(lambda v:(v*(1-v))**(rr-2)*ff(v),0,1,epsabs=1e-8,limit=200)[0]
    pred=float(answer.subs(rho,rr))*betaf(rr-1,rr-1)
    check('angular quadrature '+name,abs(val-pred)<3e-7)
eq('odd triple bulk jet averages to zero',angular_ratio(M[1]))
eq('triple endpoint odd moments are opposite',angular_ratio(M[1]*ell),-angular_ratio(M[1]*right))
eq('triple first boundary coefficient',angular_ratio(M[0]*ell),angular_ratio(M[0])/2)

# Nonconstant midpoint pair test: the endpoint counterterm cannot be omitted.
rr=-np.sqrt(2)/10; tt=2/rr; k=.3
K=2*(tt-2)
A0=np.expm1(k)/k; A1=-(1+np.exp(k))/2
def A(v):return (np.exp(k*(1-v/2))-np.exp(k*v/2))/k
def remainder(v):
    if v<1e-4:
        return sum((np.exp(k)*(-k/2)**n-(k/2)**n)*v**n/(k*float(S.factorial(n))) for n in range(2,9))
    return A(v)-A0-A1*v
target=K*(quad(lambda v:v**(rr-2)*remainder(v),0,1,epsabs=2e-10)[0]+A0/(rr-1)+A1/rr)
errors=[]
for ep in (.08,.04,.02,.01):
    raw=K*quad(lambda v:v**(rr-2)*A(v),ep,1,epsabs=1e-8)[0]
    ren=raw+K*A0*ep**(rr-1)/(rr-1)+K*A1*ep**rr/rr
    err=abs(ren-target); errors.append(err)
    print(f'DATA pair epsilon={ep:.5g} finite-part={ren:.12g} error={err:.5g}',flush=True)
check('pair including endpoint subtraction converges',all(x>y for x,y in zip(errors,errors[1:])))
check('pair remainder has predicted positive degree',abs(errors[-1]/errors[-2]-2**(-(rr+1)))<.02)

# Proper ABB source formula requires the residue at infinity of w*G.
z1,z2=S.symbols('z1 z2'); Rsum=0; Ssum=0
for zi,zj in [(z1,z2),(z2,z1)]:
    H=(a-w)**2/(w-zj)**2
    simple=S.diff(H,w).subs(w,zi)
    double=(a-zi)**2/(zi-zj)**2
    Rsum+=zi*simple; Ssum+=double
eq('two-B primitive area identity includes infinity term',Rsum+Ssum,1)

print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
