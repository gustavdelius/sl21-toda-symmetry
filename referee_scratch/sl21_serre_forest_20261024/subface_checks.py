"""Independent overlap and marginal-pole checks using explicit Wick kernels."""
import sympy as S
import numpy as np
from scipy.special import roots_jacobi
from numpy.polynomial.legendre import leggauss

passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    passed+=int(bool(ok));failed+=int(not bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)
def eq(name,a,b=0):check(name,S.expand(S.together(a-b).as_numer_denom()[0])==0)

t,r,a,R=S.symbols('t r a R',nonzero=True)
rho=2/t;eta=t-2
C,D,C2,D2,U,U2=S.symbols('C D C2 D2 U U2')
F1=-t*D-C+U;F2=-t*D2-C2+U2

# Reconstruct the triple jets directly from its three fused-tail Wick
# polynomials; the source double-pole moment gives the relative beta/t factor.
x,q=S.symbols('x q')
pos=[(1+x)/3,(1-2*x)/3,-(2-x)/3]
P2=sum(z*z for z in pos)
tail=0
for i,zi in enumerate(pos):
    js=[j for j in range(3) if j!=i];j,k=js
    zj,zk=pos[j],pos[k]
    cj=q*C+q*q*zj*C2;ck=q*C+q*q*zk*C2
    fj=cj+eta/(zj-zi)-1/(zj-zk)
    fk=ck+eta/(zk-zi)+1/(zj-zk)
    W=fj*fk+1/(zj-zk)**2
    E=1+q*F1*zi+q*q*(D2*P2+(F2+F1*F1)*zi*zi)/2
    tail-=W*E/((zj-zi)*(zk-zi))
Di=[1/S.prod((zi-zj)**2 for j,zj in enumerate(pos) if j!=i) for i,zi in enumerate(pos)]
M=[sum(zi**k*Di[i] for i,zi in enumerate(pos)) for k in range(3)]
jets=[-eta*(eta-1)*M[0],-eta*((eta-1)*F1-3*C)*M[1],
      -eta*(((eta-1)*D2/2+C2)*P2*M[0]+((eta-1)*(F2+F1*F1)/2-3*C*F1)*M[2])]
for k in range(3):
    eq(f'general-shape triple descendant from Wick tails, level {k}',S.expand(tail).coeff(q,k),jets[k])

# Embed the midpoint pair composite into the third B. The pair midpoint is
# R/3 and the remaining B is at -2R/3. There is no first odd pair correction.
pos=[S.Rational(1,3)+r/2,S.Rational(1,3)-r/2,-S.Rational(2,3)]
Di=[1/S.prod((zi-zj)**2 for j,zj in enumerate(pos) if j!=i) for i,zi in enumerate(pos)]
M=[sum(zi**k*Di[i] for i,zi in enumerate(pos)) for k in range(3)]
P2=sum(zi*zi for zi in pos)
jets=[-eta*(eta-1)*M[0],-eta*((eta-1)*F1-3*C)*M[1],
      -eta*(((eta-1)*D2/2+C2)*P2*M[0]+((eta-1)*(F2+F1*F1)/2-3*C*F1)*M[2])]
expected=[-2*eta*(eta-1),2*eta*(C-(eta-1)*F1/3),
          2*eta*(-2*C2/3+C*F1/3-(eta-1)*(D2/3+(F2+F1*F1)/18))]
for k in range(3):
    num=S.cancel(r*r*jets[k])
    eq(f'BB within BBB matches independently fused pair, jet {k}',num.subs(r,0),expected[k])
    eq(f'BB within BBB has no spurious odd midpoint pole, jet {k}',S.diff(num,r).subs(r,0))

# Independent pair embedding into the quartic scalar kernel, with the B and
# A spectators kept at arbitrary separated R and a. Strip the common vertex
# powers and r^(g-2); no g derivative is needed because it is A-independent.
pos=[r/2,-r/2,-R]
Si=sum((a-zi)**2/S.prod((zi-zj)**2 for j,zj in enumerate(pos) if j!=i) for i,zi in enumerate(pos))
logRa=rho*sum(1/(zi-a) for zi in pos)
four=S.cancel(r*r*eta/rho*(S.diff(Si,a)+logRa*Si))
leading=four.subs(r,0)
embedded=2*eta*(-a)/R*((eta-1)/(-R)+1/(-R-a))
eq('embedded pair coefficient in four-current scalar primitive',leading,embedded)
eq('four-current pair midpoint odd pole vanishes',S.diff(four,r).subs(r,0))

# UV cancellations needed when the Yukawa stays outside the current cluster.
z=[S.Rational(-3,4),S.Rational(1,5),S.Rational(9,8)]
fs=[q*C+q*q*zi*C2-sum(1/(zi-zj) for j,zj in enumerate(z) if i!=j)+1/(zi-a) for i,zi in enumerate(z)]
W=S.prod(fs)+sum(fs[k]/(z[i]-z[j])**2 for i,j,k in ((0,1,2),(0,2,1),(1,2,0)))
eq('UV BBBA scalar contraction cancels',S.expand(W).coeff(q,0))
eq('UV BBBA constant regular derivative jet cancels',S.expand(W).coeff(q,1))

# Fusion phases: each Yukawa supplies the integer-power sign and the local
# cocycle sign. The resulting composite has the phase of n original Bs.
rr=-np.sqrt(2)/7;ss=np.exp(-1j*np.pi*rr)
for n in (1,2,3):
    phaseB=-np.exp(-1j*np.pi*(n*rr-1))
    phaseA=-np.exp(-1j*np.pi*(1-n*rr))
    check(f'V{n}/B fused exchange phase',abs(phaseB-ss**n)<2e-14)
    check(f'V{n}/A fused exchange phase',abs(phaseA-ss**(-n))<2e-14)

# Compact-shape test of the only possible marginal pole. The leading scalar
# primitive P has zero A endpoints in each middle chamber. With test 1+k*a,
# integration by parts leaves -k*P; its radial factor is 1/(1+3*(g-rho)).
rr=-np.sqrt(2)/5;pp=-rr;kappa=.37
N=96
nodes,weights=roots_jacobi(N,pp,pp)
xx=(nodes+1)/2;ww=weights/(2**(2*pp+1))
dnodes,dweights=roots_jacobi(N,pp-1,pp-1)
dx=(dnodes+1)/2;dw=dweights/(2**(2*pp-1))
def primitive_shape(shape,gg,which,grid):
    zz=np.array([.5,.5-shape,-.5])
    lo,hi=(zz[2],zz[1]) if which==0 else (zz[1],zz[0])
    aa=lo+(hi-lo)*grid
    dif=zz[:,None]-aa[None,:]
    Rfac=(shape*(1-shape))**gg*np.prod(np.abs(dif)**pp,axis=0)
    Rfac=Rfac*np.exp(-1j*np.pi*rr*np.sum(dif<0,axis=0))
    DD=np.array([1/np.prod([(zz[i]-zz[j])**2 for j in range(3) if j!=i]) for i in range(3)])
    SS=np.sum((aa[None,:]-zz[:,None])**2*DD[:,None],axis=0)
    SP=2*np.sum((aa[None,:]-zz[:,None])*DD[:,None],axis=0)
    P=Rfac*SS
    derivative=Rfac*(SP+rr*np.sum(1/dif,axis=0)*SS)
    return aa,P,derivative,hi-lo
for which in (0,1):
    aa,P,derivative,length=primitive_shape(.43,rr+.17,which,dx)
    weight=(dx*(1-dx))**(pp-1)
    constant=length*np.dot(dw,derivative/weight)
    linear=length*np.dot(dw,derivative*(1+kappa*aa)/weight)
    aa,P,_,length=primitive_shape(.43,rr+.17,which,xx)
    integ=length*np.dot(ww,P/(xx*(1-xx))**pp)
    check(f'middle scalar chamber zero constant endpoint, chamber {which}',abs(constant)<3e-10,
          f'abs={abs(constant):.3g}')
    check(f'middle chamber with nonconstant test obeys Stokes, chamber {which}',abs(linear+kappa*integ)<3e-10)

xn,xw=leggauss(48);shapes=.5+.3*xn;shapeweights=.3*xw
res=[]
for delta in (.02,.01,.005,.0025,-.0025):
    gg=rr+delta;chambers=[]
    for which in (0,1):
        total=0j
        for shape,wt in zip(shapes,shapeweights):
            aa,P,_,length=primitive_shape(shape,gg,which,xx)
            total+=wt*length*np.dot(ww,P/(xx*(1-xx))**pp)
        chambers.append(-kappa*total/(1+3*delta))
    uu=np.exp(-1j*np.pi*gg);ss=np.exp(-1j*np.pi*rr)
    DD=(ss-1)*(ss+1)*(ss-uu)*(ss*ss-uu)*(uu+1)/ss**4
    defect=DD*(chambers[1]-chambers[0])
    res.append((delta,max(map(abs,chambers)),abs(defect)))
    print(f'DATA compact-shape delta={delta:.6g} max_chamber={res[-1][1]:.10g} defect={abs(defect):.10g}',flush=True)
check('compact-shape middle chambers remain bounded across diagonal',max(v[1] for v in res)<2*min(v[1] for v in res))
check('vanishing coefficient times regular nonconstant chamber tends to zero',res[3][2]<res[0][2]/5)
print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
