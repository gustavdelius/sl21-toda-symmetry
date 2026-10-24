"""Fresh Wick checks of the scalar primitive and independent-exponent extension.

No other investigation script is imported. The exponent g on BB pairs is an
auxiliary analytic regulator, not a change to the physical momentum vectors.
"""
import sympy as S
import numpy as np
from scipy.integrate import quad

passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    passed+=int(bool(ok)); failed+=int(not bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)
def eq(name,left,right=0):
    check(name,S.expand(S.together(left-right).as_numer_denom()[0])==0)

t,a,w=S.symbols('t a w',nonzero=True)
rho=2/t; beta=t-1; eta=t-2
def wick(pos,external):
    fs=[sum(c/(zi-y) for y,c in external)-sum(1/(zi-zj) for j,zj in enumerate(pos) if i!=j)
        for i,zi in enumerate(pos)]
    if len(pos)==1:return fs[0]
    if len(pos)==2:return fs[0]*fs[1]+1/(pos[0]-pos[1])**2
    return S.prod(fs)+sum(fs[k]/(pos[i]-pos[j])**2 for i,j,k in [(0,1,2),(0,2,1),(1,2,0)])

# Several asymmetric shapes; a,w,t remain symbolic. Differentiating R uses
# only its A exponent, so the identity is independent of the BB exponent g.
for n,pos in [(3,[S.Rational(-4,5),S.Rational(1,7),S.Rational(6,5)]),
              (3,[S.Rational(-2,3),S.Rational(1,3),S.Rational(7,4)]),
              (2,[S.Rational(-3,4),S.Rational(5,6)])]:
    logRa=rho*sum(1/(zi-a) for zi in pos)
    G=(a-w)**2/S.prod((w-zi)**2 for zi in pos)
    source=wick(pos,[(a,1),(w,beta)])*(a-w)/S.prod(zi-w for zi in pos)
    pref=beta*eta/rho if n==3 else beta/rho
    eq(f'Wick source A derivative n={n}, shape={pos}',source,pref*(S.diff(G,a)+logRa*G))
    for i,zi in enumerate(pos):
        rest=[zj for j,zj in enumerate(pos) if j!=i]
        Di=1/S.prod((zi-zj)**2 for zj in rest)
        Si=(a-zi)**2*Di
        tail=wick(rest,[(a,1),(zi,eta)])*(a-zi)/S.prod(zj-zi for zj in rest)
        pref=eta/rho if n==3 else 1/rho
        eq(f'Wick tail A derivative n={n}, i={i}, shape={pos}',tail,pref*(S.diff(Si,a)+logRa*Si))
    if n==3:
        eq(f'node 0 marginal scalar vanishes, shape={pos}',wick(pos,[(a,1),(w,1)]))
    if n==2:
        z1,z2=pos
        S2=sum((a-zi)**2/(z1-z2)**2 for zi in pos)
        logA=sum(1/(zi-a) for zi in pos)
        first=(2*a-z1-z2)/(z1-z2)**2
        eq('node 0 resolved ABB source cancels after integration',2*first+logA*(S2-1))
        combined=beta*(S2-1)-t*S2
        eq('ABB combined primitive includes infinity constant',combined,-S2-beta)

# Direct evaluation of all word branches at a common endpoint. Fractional
# factors are evaluated in operator-word order; no phase is supplied by hand.
rr=np.sqrt(3.); gg=np.sqrt(2.)+.5
s=np.exp(-1j*np.pi*rr); c=1+s+1/s; co=[1,-c,c,-1]
positions=[-.72,.13,.81]; levels=[1.5,.5,-.5,-1.5]
for end in (-1.,1.):
    vals=[]
    for k in range(4):
        word=['B']*(3-k)+['A']+['B']*k
        j=0; points=[]
        for i,label in enumerate(word):
            points.append((end if label=='A' else positions[j])+1e-10j*levels[i])
            j+=label=='B'
        prod=1+0j
        for i in range(4):
            for j in range(i+1,4):
                exponent=gg if word[i]==word[j]=='B' else -rr
                prod*=(points[i]-points[j])**exponent
        vals.append(prod)
    ratio=np.array(vals)/vals[0]
    expect=np.array([s**(k if end<0 else -k) for k in range(4)])
    check(f'actual word endpoint phases, end={end}',np.max(abs(ratio-expect))<1e-8)
    check(f'endpoint polynomial with independent BB exponent, end={end}',abs(np.dot(co,ratio))<3e-8)

# Equal spatial endpoints can be defined in an honest open convergence region
# g>1, rho<0 before continuing to g=rho. This is a check at fixed B positions
# of integration by parts across all three branch points, not a full integral.
rr=-np.sqrt(2.)/2; gg=np.sqrt(3.); tt=2/rr
z=np.array([-.74,.07,.65]); ss=np.exp(-1j*np.pi*rr)
co=[1,-(1+ss+1/ss),1+ss+1/ss,-1]
D=np.array([1/np.prod([(z[i]-z[j])**2 for j in range(3) if j!=i]) for i in range(3)])
Delta=complex(np.prod([z[i]-z[j] for i in range(3) for j in range(i+1,3)]))
# All B imaginary heights descend; each negative BB difference has arg +pi.
Delta_factor=np.prod([complex(z[i]-z[j])**gg for i in range(3) for j in range(i+1,3)])
def fields(aa,k):
    apos=3-k
    # Reference zi-a is above or below its cut according to the word slot.
    factors=[]
    for i,zi in enumerate(z):
        v=complex(zi-aa,0.)
        phase=0 if zi>aa else (np.pi if i<apos else -np.pi)
        factors.append(abs(v)**(-rr)*np.exp(-1j*rr*phase))
    R=ss**k*Delta_factor*np.prod(factors)
    Sch=np.dot((aa-z)**2,D)
    Sprime=2*np.dot(aa-z,D)
    derivative=R*(Sprime+rr*sum(1/(zi-aa) for zi in z)*Sch)
    return R*Sch,derivative
vals=[]
for k in range(4):
    f=lambda aa:fields(aa,k)[1]
    value=sum(quad(lambda aa:f(aa).real,lo,hi,epsabs=1e-7,limit=200)[0]
              +1j*quad(lambda aa:f(aa).imag,lo,hi,epsabs=1e-7,limit=200)[0]
              for lo,hi in zip([-1.,*z],[*z,1.]))
    boundary=fields(1.,k)[0]-fields(-1.,k)[0]
    check(f'equal-interval A integration with independent exponents, word {k}',abs(value-boundary)<3e-6,
          f'error={abs(value-boundary):.3g}')
    vals.append(value)
check('sum of integrated scalar words in convergence region',abs(np.dot(co,vals))<1e-5,
      f'absolute residual={abs(np.dot(co,vals)):.3g}')
print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
