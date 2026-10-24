"""Independent leading five-current collision reconstruction, not old-script reuse."""
import itertools
import sympy as s

t,g,rho=s.symbols('t g rho', nonzero=True)
a,w,u=s.symbols('a w u')
z=s.symbols('z0:3')
beta,eta=t-1,t-2

def ck(name, x):
    v=s.factor(s.together(x))
    assert v == 0, (name,v)
    print('PASS '+name)

# Normal-ordered scalar coefficient from all Wick matchings of three C's.
linear=[1/(zi-a)+beta/(zi-w)-sum(1/(zi-zj) for j,zj in enumerate(z) if j!=i) for i,zi in enumerate(z)]
wick=s.prod(linear)+sum(linear[k]/(z[i]-z[j])**2 for i,j,k in [(0,1,2),(0,2,1),(1,2,0)])
# Scalar vertex factor = R * (a-w)^2/prod(zi-w), with R carrying
# like-current powers and A/B powers; normal fields not silently differentiated.
G=(a-w)**2/s.prod((w-zi)**2 for zi in z)
vertex_extra=(a-w)/s.prod(zi-w for zi in z)
logRa=-2/t*sum(1/(a-zi) for zi in z)
candidate=beta*eta/(2/t)*(s.diff(G,a)+logRa*G)
ck('three B Gaussian Wick polynomial is minority derivative',vertex_extra*wick-candidate)

# With barB at u, its scalar Y2 contraction is beta/(u-wbar)^2.
wb=s.symbols('wb')
ck('barB scalar factor is an independent derivative',s.diff(-beta/(u-wb),u)-beta/(u-wb)**2)
# Since R and G contain no u, the scalar five-current/source coefficient
# is a double derivative. This alone says nothing about its proper faces.
ck('five-current scalar double derivative coefficient',
   candidate*beta/(u-wb)**2-s.diff(-beta*candidate/(u-wb),u))

# Independent even-root minority derivative, plus barA second-order pole.
zz=s.symbols('zz')
aa=s.symbols('aa0:3')
G_A=(zz-w)**2/s.prod((w-ai)**2 for ai in aa)
logRz=-2/t*sum(1/(zz-ai) for ai in aa)
L_B=sum(1/(zz-ai) for ai in aa)-t/(zz-w)
ck('even-root holomorphic minority derivative',G_A*L_B+(s.diff(G_A,zz)+logRz*G_A)/(2/t))
ck('barA even-root pole is a derivative',s.diff(-1/(u-wb),u)-1/(u-wb)**2)

# All source-containing subsets of B B B A barB; no current-current
# hol/antihol power, and derivative contractions counted from the actual set.
labels=['B','B','B','A','D']
rows=[]
for n in range(1,6):
    for subset in itertools.combinations(range(5),n):
        types=[labels[i] for i in subset]
        nb,na,nd=types.count('B'),types.count('A'),types.count('D')
        pair=g*nb*(nb-1)/2-rho*nb*na
        source=na-nb-nd
        derivatives=nb+nd
        measure=n+1
        degree=s.expand(pair+source-derivatives+measure)
        rows.append((tuple(types),degree))
ck('BBbarB source degree', next(d for x,d in rows if x==('B','B','D'))-(g-2))
ck('five-current source degree',rows[-1][1]-(3*(g-rho)-1))
ck('double-primitive degree gain',rows[-1][1]+2-(3*(g-rho)+1))
for types,degree in sorted(set(rows),key=lambda x:(len(x[0]),x[0])):
    print('INFO subset',types,'relative degree',degree)

# A generic analytic regulator is not multiplicative under scalar finite part.
# This witnesses why n-ary uniqueness cannot replace a composition proof.
e=s.symbols('e')
fp=lambda expr:s.expand(expr).coeff(e,0)
assert fp(e)*fp(1/e)==0 and fp(e/e)==1
print('PASS finite-part operation is not an algebra homomorphism (explicit counterexample)')
print('INFO Scalar double derivatives require extension through mixed/like proper faces; no associativity conclusion here.')
