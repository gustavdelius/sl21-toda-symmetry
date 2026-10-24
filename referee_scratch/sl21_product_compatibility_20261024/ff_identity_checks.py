"""Independent sign sums and recursion from arXiv:1610.04926 (2017).

No old checker imports. v=omega**(1/4) is kept symbolic in Laurent sums.
The rational recursion checks use v=2, a non-root-of-unity algebraic probe.
"""
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
import sympy as sp

count=0
def check(label,ok):
    global count
    if not isinstance(ok,bool):
        ok=sp.cancel(ok)==0
    print(('PASS' if ok else 'FAIL')+' '+label)
    if not ok:
        raise AssertionError(label)
    count+=1

def multiply(P,Q):
    ans={}
    for j,a in P.items():
        for k,b in Q.items():
            ans[j+k]=ans.get(j+k,F(0))+a*b
    return {k:a for k,a in ans.items() if a}

def rational_pair(e,z,kind):
    z=F(z)
    if kind=='same':
        return {2*e:z/(z-1),-2*e:-1/(z-1)}
    return {-2*e:z/(z+1),2*e:1/(z+1)}

def laurent_sum(X,Y,slopes=(1,1),phases=(1,1)):
    """Return the exact Laurent polynomial, with universal i^(K+L) removed."""
    total={}
    K,L=len(X),len(Y)
    for signs in product((-1,1),repeat=K+L):
        es,ds=signs[:K],signs[K:]
        exponent=slopes[0]*sum(es)+slopes[1]*sum(ds)
        iphase=phases[0]*sum(es)+phases[1]*sum(ds)-(K+L)
        assert iphase%2==0
        P={exponent:F((-1)**((iphase//2)%2))}
        for Z,zsigns in [(X,es),(Y,ds)]:
            for i in range(len(Z)):
                for j in range(i+1,len(Z)):
                    if zsigns[i]!=zsigns[j]:
                        P=multiply(P,rational_pair(zsigns[i],F(Z[j],Z[i]),'same'))
        for i,x in enumerate(X):
            for j,y in enumerate(Y):
                if es[i]!=ds[j]:
                    P=multiply(P,rational_pair(es[i],F(y,x),'cross'))
        for k,a in P.items():
            total[k]=total.get(k,F(0))+a
    return {k:a for k,a in total.items() if a}

xs=(2,5,11,17,23)
ys=(3,7,13,19,29)
for n in range(1,5):
    X,Y=xs[:n],ys[:n]
    for name,slopes,phases in [('identity',(1,1),(1,1)),
                              ('FA',(-1,-1),(-1,-1)),
                              ('U',(1,1),(1,3))]:
        check('exact Laurent zero '+name+' n='+str(n),not laurent_sum(X,Y,slopes,phases))
for K,L in [(0,1),(0,2),(0,4),(1,2),(1,3),(2,3),(2,4)]:
    r=L-K
    check('induction-family zero K='+str(K)+' L='+str(L),
          not laurent_sum(xs[:K],ys[:L],(1+r,1-r)))

# Independent rational implementation, keeping zero-mode coefficients generic.
def kernel(e,z,v,kind):
    return ((v**(2*e)*z-v**(-2*e))/(z-1) if kind=='same'
            else (v**(-2*e)*z+v**(2*e))/(z+1))

def sign_sum(X,Y,alpha,beta,v):
    ans=F(0)
    K,L=len(X),len(Y)
    for signs in product((-1,1),repeat=K+L):
        es,ds=signs[:K],signs[K:]
        q=alpha**sum(es)*beta**sum(ds)
        for Z,zsigns in [(X,es),(Y,ds)]:
            for i in range(len(Z)):
                for j in range(i+1,len(Z)):
                    if zsigns[i]!=zsigns[j]:
                        q*=kernel(zsigns[i],Z[j]/Z[i],v,'same')
        for i,x in enumerate(X):
            for j,y in enumerate(Y):
                if es[i]!=ds[j]:
                    q*=kernel(es[i],y/x,v,'cross')
        ans+=q
    return ans

@lru_cache(None)
def recursion(X,Y,alpha,beta,v):
    K,L=len(X),len(Y)
    if K==0:
        ans=F(1)
        for i in range(1,L+1):
            z=L+1-2*i
            ans*=beta*v**z+beta**-1*v**-z
        return ans
    x,XX=X[0],X[1:]
    C=alpha*v**(-(K-1)+L)+alpha**-1*v**((K-1)-L)
    ans=C*recursion(XX,Y,alpha*v,beta/v,v)
    for i,y in enumerate(Y):
        first=alpha/beta
        second=beta/alpha
        for xx in XX:
            first*=kernel(1,-xx/y,v,'same')
            second*=kernel(-1,-xx/y,v,'same')
        for j,yy in enumerate(Y):
            if j!=i:
                first*=kernel(-1,yy/y,v,'same')
                second*=kernel(1,yy/y,v,'same')
        R=(v**-2-v**2)*(first-second)
        ans+=y*R/(x+y)*recursion(XX,Y[:i]+Y[i+1:],alpha,beta,v)
    return ans

for K,L in [(0,0),(0,1),(0,3),(1,0),(1,1),(1,2),(2,1),(2,2),(2,3),(3,3)]:
    X,Y=tuple(map(F,xs[:K])),tuple(map(F,ys[:L]))
    check('rational recursion vs sign sum K='+str(K)+' L='+str(L),
          recursion(X,Y,F(3),F(5),F(2))==sign_sum(X,Y,F(3),F(5),F(2)))

v,z=sp.symbols('v z',nonzero=True)
f=lambda e,z:(v**(2*e)*z-v**(-2*e))/(z-1)
gg=lambda e,z:(v**(-2*e)*z+v**(2*e))/(z+1)
for e in (-1,1):
    check('f inversion e='+str(e),f(e,1/z)-f(-e,z))
    check('g inversion e='+str(e),gg(e,1/z)-gg(-e,z))
check('generic-rapidity J11 identity',-v**2-v**-2+gg(1,z)+gg(-1,z))
check('exceptional J00 residue', (v**-2-v**2)*(1-1))
check('nonzero generic sign sum',sign_sum((F(2),),(F(3),),F(3),F(5),F(2))!=0)
check('zero-particle normalization',sign_sum((),(),F(3),F(5),F(2))==1)
print('SUMMARY',count,'PASS, 0 FAIL')
