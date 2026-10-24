"""Independent first-order shuffle, tail-cocycle, and forest compatibility audit.

No previous check is imported. Fractional powers use actual operator order.
Nonconstant coherent backgrounds test the complete Gaussian Wick polynomial.
"""
from itertools import permutations
import sympy as S
import numpy as np

passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    passed+=int(bool(ok));failed+=int(not bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)
def eq(name,a,b=0):check(name,S.expand(S.together(a-b).as_numer_denom()[0])==0)

u,s=S.symbols('u s',nonzero=True)
words=['BBBA','BBAB','BABB','ABBB']
M=S.zeros(4)
allperms=list(permutations(range(4)))
for j,word in enumerate(words):
    for perm in allperms:
        rank={label:k for k,label in enumerate(perm)}
        weight=S.Integer(1)
        for i in range(4):
            for k in range(i+1,4):
                if rank[i]>rank[k]:weight*=u if word[i]==word[k] else 1/s
        M[words.index(''.join(word[i] for i in perm)),j]+=weight
c=S.Matrix([1,-1-s-1/s,1+s+1/s,-1])
defect=(s-1)*(s+1)*(s-u)*(s*s-u)*(u+1)/s**4
for row in range(4):eq(f'two-exponent permutation defect row {row}',(M*c)[row],[0,-defect,defect,0][row])
eq('generic diagonal null vector',sum(v*v for v in (M*c).subs(u,s)))
print('DATA defect D =',S.factor(defect),flush=True)

# Apply the same chamber-linear functional to every local Taylor coefficient.
# Symbols stand for independently varying background and endpoint jets.
for face,levels in [('BB',2),('BBB',3),('AB',1),('ABB',2)]:
    for level in range(levels):
        jets=S.Matrix(1,4,S.symbols(f'{face}_{level}_0:4'))
        result=(jets*M*c)[0]
        eq(f'embedded {face} subtraction with arbitrary chamber jet {level}',
           result,defect*(jets[2]-jets[1]))
        eq(f'embedded {face} subtraction vanishes on diagonal, jet {level}',result.subs(u,s))

def wick(values,points):
    n=len(values)
    if n==0:return 1.+0j
    if n==1:return values[0]
    if n==2:return values[0]*values[1]+1/(points[0]-points[1])**2
    return np.prod(values)+sum(values[k]/(points[i]-points[j])**2 for i,j,k in ((0,1,2),(0,2,1),(1,2,0)))

def amplitude(word,coords,g,rho,node,channel):
    t=2/rho; beta=1 if node==0 else t-1
    source=.21+.73j
    q=np.array([1 if label=='B' else -1 for label in word])
    # C(x) and d.X(x) vary independently; a2.X=-t*d.X-a0.X.
    C=lambda x:.3+.2*x+.11*x*x
    phi=lambda x:.19*x+.037*x*x
    phi0=lambda x:.3*x+.1*x*x+.11*x**3/3
    antihol=lambda x:(.13*x+.04*x*x) if node==2 else (.07*x-.02*x*x)
    screen=lambda x:np.exp((phi0(x) if node==0 else -t*phi(x)-phi0(x))+antihol(np.conj(x)))
    base=np.exp(sum(q[i]*phi(x) for i,x in enumerate(coords)))
    tails=[None] if channel=='action' else [i for i,label in enumerate(word) if label=='B']
    total=0j
    for tail in tails:
        derivatives=[i for i,label in enumerate(word) if label=='B' and i!=tail]
        vertex=1.+0j
        for i in range(4):
            for j in range(i+1,4):
                exponent=g if word[i]==word[j] else -rho
                if tail in (i,j):exponent-=q[j if tail==i else i]
                vertex*=complex(coords[i]-coords[j])**exponent
        vals=[]
        for i in derivatives:
            v=C(coords[i])
            for j in range(4):
                if j==i:continue
                v+=(t-2 if j==tail else -q[j])/(coords[i]-coords[j])
            if channel=='action':v+=beta/(coords[i]-source)
            vals.append(v)
        W=wick(vals,[coords[i] for i in derivatives])
        if channel=='action':
            # Net current/source cocycle is +1 for either quartic.
            vertex*=np.prod([(coords[i]-source)**(-q[i]) for i in range(4)])
            total+=-base*vertex*W*screen(source)
        else:
            total+=-t*(-1)**tail*base*vertex*W*screen(coords[tail])
    return total

for setno,(gg,rr) in enumerate([(np.sqrt(2.),-np.sqrt(3.)/5),
                              (np.sqrt(3.)+.1j,-np.sqrt(2.)/8+.07j)]):
    uu=np.exp(-1j*np.pi*gg);ss=np.exp(-1j*np.pi*rr)
    coef=np.array([1,-1-ss-1/ss,1+ss+1/ss,-1])
    for major in ('B','A'):
        wordset=words if major=='B' else [w.translate(str.maketrans('AB','BA')) for w in words]
        xs=np.array([-.91,-.23,.42,.88])
        for node,channel in ((0,'action'),(2,'action'),(2,'tail')):
            worst=0.; totals=np.zeros(4,complex)
            canonical=np.array([amplitude(w,xs,gg,rr,node,channel) for w in wordset])
            for col,word in enumerate(wordset):
                for perm in allperms:
                    rank={label:k for k,label in enumerate(perm)}
                    coords=np.array([xs[rank[i]] for i in range(4)])
                    row=wordset.index(''.join(word[i] for i in perm))
                    factor=1.+0j
                    for i in range(4):
                        for j in range(i+1,4):
                            if rank[i]>rank[j]:factor*=uu if word[i]==word[j] else 1/ss
                    actual=amplitude(word,coords,gg,rr,node,channel)
                    expected=factor*canonical[row]
                    worst=max(worst,abs(actual-expected)/(1+abs(expected)))
                    totals[col]+=actual
            check(f'complete Wick/cocycle covariance {major}^3, node {node}, {channel}, set {setno}',worst<2e-12,
                  f'max scaled error={worst:.3g}')
            DD=(ss-1)*(ss+1)*(ss-uu)*(ss*ss-uu)*(uu+1)/ss**4
            expected=DD*(canonical[2]-canonical[1])
            error=abs(np.dot(coef,totals)-expected)/(1+abs(expected))
            check(f'full nonconstant-background word sum {major}^3, node {node}, {channel}, set {setno}',error<5e-11,
                  f'scaled error={error:.3g}')

# Audit the potentially marginal lowest source coefficient from its Wick
# definition, with unrestricted a,w,t and asymmetric rational B positions.
a,w,t=S.symbols('a w t',nonzero=True);rho=2/t;beta=t-1;eta=t-2
z=[S.Rational(-5,6),S.Rational(1,8),S.Rational(7,5)]
fi=[-sum(1/(zi-zj) for j,zj in enumerate(z) if i!=j)+1/(zi-a)+beta/(zi-w) for i,zi in enumerate(z)]
W=S.prod(fi)+sum(fi[k]/(z[i]-z[j])**2 for i,j,k in ((0,1,2),(0,2,1),(1,2,0)))
G=(a-w)**2/S.prod((w-zi)**2 for zi in z)
logR=rho*sum(1/(zi-a) for zi in z)
eq('independent scalar source exact derivative',W*(a-w)/S.prod(zi-w for zi in z),beta*eta/rho*(S.diff(G,a)+logR*G))
eq('node 0 primitive scalar zero',W.subs(t,2))

print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
