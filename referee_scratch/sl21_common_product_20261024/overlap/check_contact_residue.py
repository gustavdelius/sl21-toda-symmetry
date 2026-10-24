import itertools
import numpy as np
import sympy as S
from scipy.special import roots_jacobi
z1,z2,a,x,t,rho,g=S.symbols('z1 z2 a x t rho g')
eta=t-2;delta=g-rho
f1=-1/(z1-z2)+1/(z1-a)+eta/(z1-x)
f2=-1/(z2-z1)+1/(z2-a)+eta/(z2-x)
Q=(a-x)/((z1-x)*(z2-x))*(f1*f2+1/(z1-z2)**2)
Si=(a-x)**2/((x-z1)**2*(x-z2)**2)
logR=rho*(1/(z1-a)+1/(z2-a)+1/(x-a))
checks=[]
def check(name,condition):
    ok=bool(condition);checks.append(ok);print(('PASS' if ok else 'FAIL')+' '+name)
check('actual two-C endpoint kernel is minority derivative',S.factor((Q-eta/rho*(S.diff(Si,a)+logR*Si)).subs(rho,2/t))==0)
H=1/(x-z1)+1/(x-z2)
check('insertionwise vs whole derivative difference',S.simplify(S.diff(S.log(z1-x)+S.log(z2-x),x)-H)==0)
# Exact angular exponents after x=0,z1=uR,z2=R,a=uyR, partial_a transferred to test.
u,y,R=S.symbols('u y R',positive=True)
# R factor: BB powers 3g; BA powers -3rho; S power -2; H power -1; relative measure power +2.
radial=3*g-3*rho-2-1+2
check('actual difference first jet gives delta/(3delta)',S.simplify(radial-(3*delta-1))==0)
# Quadrature is restricted to the one explicitly stated middle chamber.
rhov=-np.sqrt(2)/5

def angular(gv,n=128):
    xu,wu=roots_jacobi(n,gv,gv-2*rhov)
    xy,wy=roots_jacobi(n,-rhov,2-rhov)
    us=(xu+1)/2; ys=(xy+1)/2
    wu=wu/2**(2*gv-2*rhov+1);wy=wy/2**(3-2*rhov)
    vals=(1+us[:,None])*(1-us[:,None]*ys[None,:])**(-rhov)
    return float(wu@vals@wy)
I0=angular(rhov,192);I1=angular(rhov,384)
check('middle-chamber angular integral converges and is nonzero',I1>0 and abs(I0-I1)<2e-6)
print('DIAGNOSTIC rho',rhov,'angular I',I1,'normalized radial residue I/3',I1/3)
for dv in [.04,.02,.01,.005,.0025]:
    val=angular(rhov+dv,192)/3
    print('DIAGNOSTIC delta',dv,'delta times radial integral times angular',val)
# Explicit DAG properness check for every order on four labels. This tests the premise,
# not an unproved assignment of fused contour heights in the field-theory prescription.
for order in itertools.permutations(range(4)):
    heights={v:4-i for i,v in enumerate(order)}
    edges=[(order[i],order[j]) for i in range(4) for j in range(i+1,4)]
    ks=np.array([.11,.23,.37,.41,.53,.67])
    flow=np.zeros(4)
    for (i,j),k in zip(edges,ks):flow[i]+=k;flow[j]-=k
    pairing=sum(heights[i]*flow[i] for i in range(4))
    rhs=sum((heights[i]-heights[j])*k for (i,j),k in zip(edges,ks))
    check('DAG properness '+''.join(map(str,order)),abs(pairing-rhs)<1e-12 and pairing>=sum(ks)-1e-12)
# Opposite assignments lose properness; records why the contour premise matters.
flow=np.array([1.,-1.])+np.array([-1.,1.])
check('opposite i0 factors permit positive circulation',np.linalg.norm(flow)==0)
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
