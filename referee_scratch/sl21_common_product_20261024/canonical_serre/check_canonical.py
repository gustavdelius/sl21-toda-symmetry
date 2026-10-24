import cmath, itertools, math
from functools import lru_cache
import sympy as S
from pathlib import Path
checks=[]
def test(name,condition):
    ok=bool(condition);checks.append(ok);print(('PASS' if ok else 'FAIL')+' '+name)
# Fresh Wick records: q,m0,m2,derivativeC,original current label.
def rec(label,kind='bare'):
    if label=='A':return (-1,-1,-1,False,label) if kind=='F' else (-1,0,0,False,label)
    return {'bare':(1,0,0,True,label),'K':(1,0,1,False,label),'L':(1,1,1,True,label)}[kind]
def momenta(v,t):
    q,m,n,_,_=v
    return (q,m,n)
def power(x,p,anti=False):
    return cmath.exp(p*(math.log(abs(x))+(1j*(-math.pi if anti else math.pi) if x<0 else 0)))
def wick(vertices,x,t,g):
    rho=2/t; val=1+0j
    for i in range(len(x)):
        qi,mi,ni,di,li=vertices[i]
        # Nonconstant coherent normal fields; a2=-t*d-a0 separately in both chiralities.
        D=lambda y:.071*y+.013*y*y+.004*y*y*y
        P=lambda y:.11*y+.019*y*y+.002*y*y*y
        Db=lambda y:-.033*y+.007*y*y
        Pb=lambda y:.043*y-.011*y*y
        val*=cmath.exp(qi*D(x[i])+mi*P(x[i])+ni*(-t*D(x[i])-P(x[i]))+mi*Pb(x[i])+ni*(-t*Db(x[i])-Pb(x[i])))
        for j in range(i+1,len(x)):
            qj,mj,nj,dj,lj=vertices[j]
            anti=mi*mj+ni*nj+(t-1)*(mi*nj+ni*mj)
            hol=rho*qi*qj-qi*(mj+nj)-qj*(mi+ni)+anti+(g-rho)*int(li==lj)
            val*=(-1)**(qi*(mj+nj))*power(x[i]-x[j],hol)*power(x[i]-x[j],anti,True)
    ds=tuple(i for i,v in enumerate(vertices) if v[3])
    linear={i:.11+.038*x[i]+.006*x[i]*x[i]+sum((-v[0]+v[1]+(t-1)*v[2])/(x[i]-x[j]) for j,v in enumerate(vertices) if j!=i) for i in ds}
    @lru_cache(None)
    def poly(ids):
        if not ids:return 1+0j
        i,*rest=ids;ans=linear[i]*poly(tuple(rest))
        for j in rest:ans+=poly(tuple(k for k in rest if k!=j))/(x[i]-x[j])**2
        return ans
    return val*poly(ds)
def channels(word):
    yield 'zero',[rec(c) for c in word]
    bi=[i for i,c in enumerate(word) if c=='B'];ai=[i for i,c in enumerate(word) if c=='A']
    for i in bi:
        vs=[rec(c) for c in word];vs[i]=rec('B','K');yield 'oneK'+str(i),vs
    for ij in itertools.combinations(bi,2):
        vs=[rec(c) for c in word]
        for i in ij:vs[i]=rec('B','K')
        yield 'twoK'+str(ij),vs
    for i in ai:
        vs=[rec(c) for c in word];vs[i]=rec('A','F');yield 'oneF'+str(i),vs
    for i in bi:
        vs=[rec(c) for c in word];vs[i]=rec('B','L');yield 'oneL'+str(i),vs
# Actual phase comparison for every permutation and tail position, full Wick polynomial retained.
xs=[-.9,-.17,.63,1.41]
for t,g in [(-math.sqrt(2),6.71),(math.sqrt(5),1.17),(math.sqrt(5),2/math.sqrt(5))]:
    rho=2/t;u=cmath.exp(-1j*math.pi*g);s=cmath.exp(-1j*math.pi*rho)
    for majority,minority in [('B','A'),('A','B')]:
        for pos in range(4):
            word=majority*pos+minority+majority*(3-pos)
            for cname,vs in channels(word):
                err=0
                for perm in itertools.permutations(range(4)):
                    x=[xs[perm[i]] for i in range(4)]
                    order=sorted(range(4),key=lambda i:x[i])
                    base=wick([vs[i] for i in order],sorted(x),t,g)
                    phase=1+0j
                    for i in range(4):
                        for j in range(i+1,4):
                            if perm[i]>perm[j]:phase*=u if word[i]==word[j] else 1/s
                    actual=wick(vs,x,t,g)
                    err=max(err,abs(actual-phase*base)/max(1,abs(actual),abs(base)))
                test('full Wick phase '+word+'/'+cname+' t='+str(round(t,5))+' g='+str(round(g,5)),err<3e-11)
# Exact labelled shuffle, independent enumeration.
u,s=S.symbols('u s',nonzero=True)
words=['BBBA','BBAB','BABB','ABBB'];M=S.zeros(4)
for col,word in enumerate(words):
    for perm in itertools.permutations(range(4)):
        spatial=''.join(word[i] for i in sorted(range(4),key=lambda i:perm[i]))
        phase=S.Integer(1)
        for i in range(4):
            for j in range(i+1,4):
                if perm[i]>perm[j]:phase*=u if word[i]==word[j] else 1/s
        M[words.index(spatial),col]+=phase
c=1+s+1/s;vec=S.Matrix([1,-c,c,-1]);out=M*vec
D=(s-1)*(s+1)*(s-u)*(s*s-u)*(u+1)/s**4
for i,w in enumerate(words):test('exact shuffle coefficient '+w,S.factor(out[i]-[0,-D,D,0][i])==0)
# Generic source-free and individual-tail primitive identities at asymmetric positions.
t,rho,a,z=S.symbols('t rho a z',nonzero=True)
eta=t-2; gamma=t/(1-t)
for pts in [[S.Rational(-2),S.Rational(1,3),S.Rational(3)], [S.Rational(-5,3),S.Rational(1,2),S.Rational(11,4)]]:
    # Bare B^3 A: fully contracted scalar and one constant C insertion vanish.
    f=[-sum(1/(pts[i]-pts[j]) for j in range(3) if j!=i)+1/(pts[i]-a) for i in range(3)]
    scalar=S.prod(f)+sum(f[3-i-j]/(pts[i]-pts[j])**2 for i in range(3) for j in range(i+1,3))
    oneC=sum(f[i]*f[j]+1/(pts[i]-pts[j])**2 for i in range(3) for j in range(i+1,3))
    test('bare B3A scalar cancellation '+str(pts),S.factor(scalar)==0)
    test('bare B3A constant-C cancellation '+str(pts),S.factor(oneC)==0)
    logRa=rho*sum(1/(v-a) for v in pts)
    for i in range(3):
        others=[j for j in range(3) if j!=i];j,k=others
        fi=-1/(pts[j]-pts[k])+1/(pts[j]-a)+eta/(pts[j]-pts[i]);fk=-1/(pts[k]-pts[j])+1/(pts[k]-a)+eta/(pts[k]-pts[i])
        extra=(a-pts[i])/((pts[j]-pts[i])*(pts[k]-pts[i]))
        Q=extra*(fi*fk+1/(pts[j]-pts[k])**2)
        Si=(a-pts[i])**2/S.prod((pts[i]-pts[jj])**2 for jj in others)
        test('oneK scalar individual derivative '+str(i)+' '+str(pts),S.factor((Q-eta/rho*(S.diff(Si,a)+logRa*Si)).subs(rho,2/t))==0)
        # AAAB, one FA: minority B derivative with one altered A momentum.
        SiA=(z-pts[i])**2/S.prod((pts[i]-pts[jj])**2 for jj in others)
        fB=sum(1/(z-v) for v in pts)-t/(z-pts[i])
        logRz=-rho*sum(1/(z-v) for v in pts)
        test('oneF scalar individual derivative '+str(i)+' '+str(pts),S.factor((SiA*fB+(S.diff(SiA,z)+logRz*SiA)/rho).subs(rho,2/t))==0)
# Degress generated from actual records and contracted derivatives, all relevant subsets.
g,rho,t=S.symbols('g rho t')
def degree(vs):
    e=0
    for i,v in enumerate(vs):
        qi,mi,ni,di,li=v
        for qj,mj,nj,dj,lj in vs[i+1:]:
            anti=mi*mj+ni*nj+(t-1)*(mi*nj+ni*mj)
            e+=rho*qi*qj-qi*(mj+nj)-qj*(mi+ni)+2*anti+(g-rho)*int(li==lj)
    return S.expand(e-sum(v[3] for v in vs)+len(vs)-1)
expected={'B':{'oneK':0,'twoK':2,'oneF':6,'oneL':-2},'A':{'oneK':6,'oneF':0,'oneL':8}}
for maj,minr in [('B','A'),('A','B')]:
    word=maj*3+minr
    for name,vs in channels(word):
        kind=next((k for k in expected[maj] if name.startswith(k)),None)
        if kind:test('full degree '+word+'/'+name,S.simplify(degree(vs)-(3*(g-rho)+expected[maj][kind]))==0)
# Proper subfaces at rho in (2/3,1), check endpoints of linear rho powers after physical g=rho.
for maj,minr in [('B','A'),('A','B')]:
    word=maj*3+minr
    for name,vs in channels(word):
        if name=='zero' or name.startswith('oneL'):continue
        degrees=set()
        for count in [2,3]:
            for ids in itertools.combinations(range(4),count):
                sub=[vs[i] for i in ids]
                if all(v[1:3]==(0,0) for v in sub):continue
                physical=S.factor(degree(sub).subs(g,rho))
                degrees.add(str(physical))
                # A generic pole may be inherited, but no proper face is identically marginal.
                test('proper degree not identically zero '+word+'/'+name+'/'+str(ids),physical!=0)
        print('DIAGNOSTIC proper degrees',word,name,sorted(degrees))
# The actual cutoff improvement vanishes uniformly on the physical open patch 2<t<3.
x=S.symbols('x',real=True)
for maj,minr in [('B','A'),('A','B')]:
    word=maj*3+minr
    for name,vs in channels(word):
        if not name.startswith('oneL'):continue
        li=next(i for i,v in enumerate(vs) if v[1:3]==(1,1))
        for count in [2,3,4]:
            for ids in itertools.combinations(range(4),count):
                if li not in ids:continue
                dg=degree([vs[i] for i in ids])+2*t
                polynomial=S.Poly(S.expand(S.cancel(t*dg.subs(g,rho).subs(rho,2/t)).subs(t,x+2)),x)
                # Positive monomial coefficients on x in [0,1] suffice here.
                test('cutoff L face vanishes '+word+'/'+name+'/'+str(ids),all(v>=0 for v in polynomial.all_coeffs()) and polynomial.eval(0)>0)
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
