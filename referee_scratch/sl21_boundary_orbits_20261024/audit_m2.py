"""Independent root/Wick, Ward and new mass-two collision checks.

No earlier script is imported. Exact inequalities concern a formal complex
continuation near 2<t<3, not a physical real-coupling domain.
"""
from itertools import combinations
import sympy as S

checks = 0
def check(name, value):
    global checks
    if not bool(value):
        print('FAIL', name)
        raise AssertionError(name)
    checks += 1
    print('PASS', name)

def zero(e):
    return S.simplify(e) == 0

t, rho, g, h = S.symbols('t rho g h', nonzero=True)
a0 = S.Matrix([S.I*S.sqrt(t/2), S.sqrt(t/2), S.I])
a2 = S.Matrix([S.I*S.sqrt(t/2), -S.sqrt(t/2), -S.I])
T = a0+a2
d = -T/t
check('screening Gram matrix', a0.dot(a0)==-1 and a2.dot(a2)==-1
      and zero(a0.dot(a2)-(1-t)))
check('pair-dual products', zero(d.dot(d)+2/t) and zero(d.dot(a0)-1)
      and zero(d.dot(a2)-1))
check('screening-first node-zero residue vanishes', a0+a0.dot(a0)*a0 == S.zeros(3,1))
check('screening-first node-two residue is -t(d+a2)',
      all(zero(x) for x in a0+a0.dot(a2)*a2+t*(d+a2)))

def wick(labels, pair):
    if not labels: return S.Integer(1)
    return sum((-1)**(j-1)*pair(labels[0],labels[j])*
               wick(labels[1:j]+labels[j+1:],pair)
               for j in range(1,len(labels)))

def fermion_pair(a,b):
    # Numeric labels: chirality, species, point. Only distinct points contract.
    return int(a[0]==b[0] and a[1]!=b[1] and a[2]!=b[2])
YY=[('L',0,0),('R',0,0),('L',2,1),('R',2,1)]
jj=[('L',0,0),('L',2,0),('L',0,1),('L',2,1)]
check('full Yukawa fusion sign includes i squared and one crossing',
      S.I**2*wick(YY,fermion_pair)==1)
check('fermion current level agrees with i partial h',wick(jj,fermion_pair)==1)

C,D,C2,D2 = S.symbols('C D C2 D2')
N=(C2+C*C-(-t*D2-C2+(-t*D-C)**2))/2
simple=C2+C*D+t*(D2+D*D)/2
diff_double=C2+t*D2+(C+t*D)*(1-t)*D
diff2_triple=t*((1-t)*D2+(1-t)**2*D*D)
# The action expansion contributes -kappa/pi. The pole combination alone has
# the opposite sign from the Ward response multiplying +kappa.
RP=-S.expand(simple-diff_double+diff2_triple/2)
check('all three potential poles including the minus action sign',
      zero(RP-(-t*C*D+t*t*D2/2-t**3*D*D/2)))
check('double-Yukawa contact minus potential gives partial L',
      zero(N-RP/t-(C2+(1-t)*C*D)))

eps,R = S.symbols('eps R', positive=True)
r = S.symbols('r', positive=True)
radial=S.pi*(R**(2*h)-eps**(2*h))/h
check('unlike-source radial normalization',zero(S.diff(radial,R)-2*S.pi*R**(2*h-1)))
check('source lower disk coefficient gives kappa_h',
      zero(radial/S.pi**2+eps**(2*h)/(S.pi*h)-R**(2*h)/(S.pi*h)))
check('tail lower disk coefficient is t/h, not one',
      zero(t*radial/S.pi+t*eps**(2*h)/h-t*R**(2*h)/h))

# Reconstruct every source-containing subset degree from vertex contractions,
# integration dimensions, derivative count and explicit cutoff powers.
# Token = (q,m0,m2,has_C,is_current_slot).
A=(-1,0,0,0,1); B=(1,0,0,1,1)
Y0=(0,1,0,0,0); Y2=(0,0,1,0,0)
K=(1,0,1,0,1); L=(1,1,1,1,1)
P=(0,1,1,0,0); V=(0,-1,-1,0,0); FA=(-1,-1,-1,0,1)
def pair_power(x,y):
    q,m,n,_,slot=x; Q,M,N,_,Slot=y
    screen=m*M+n*N+(h-1)*(m*N+n*M)
    return rho*q*Q+(g-rho)*int(slot and Slot and q==Q)-q*(M+N)-Q*(m+n)+2*screen

def degree(tokens,cutoff=0):
    nc=sum(x[4] for x in tokens)
    ns=len(tokens)-nc
    measure=nc+2*ns-(1 if nc else 2)
    return S.expand(sum(pair_power(x,y) for x,y in combinations(tokens,2))
                    -sum(x[3] for x in tokens)+measure+cutoff)

def expected(a,b,channel):
    E=g*(a*(a-1)+b*(b-1))/2-rho*a*b
    if channel=='unlike': return E+3*a-2*b+2*h+1
    if channel=='equal': return E+3*a-2*b+5
    return E-a+2*b+1

def positive_patch(expr):
    # Exact quadratic Bernstein representation on 2<=t<=3.
    x=S.symbols('x',real=True)
    poly=S.Poly(S.cancel(t*expr.subs({g:2/t,rho:2/t,h:t})).subs(t,2+x),x)
    n=poly.degree()
    bern=[S.simplify(sum(poly.nth(k)*S.binomial(j,k)/S.binomial(n,k)
                        for k in range(j+1))) for j in range(n+1)]
    return all(v>=0 for v in bern) and any(v>0 for v in bern)

subsets=sorted(set((a,b) for na,nb in [(1,3),(3,1)]
                   for a in range(na+1) for b in range(nb+1) if a+b))
for a,b in subsets:
    cur=[A]*a+[B]*b
    for name,sources in [('unlike',[Y0,Y2]),('equal',[Y0,Y0]),('even',[V])]:
        deg=degree(cur+sources)
        check(f'raw subset a={a} b={b} {name}',zero(deg-expected(a,b,name)))
        if name!='even':
            check(f'positive physical-patch degree a={a} b={b} {name}',positive_patch(deg))
    if b:
        other=[A]*a+[B]*(b-1)
        check(f'fused K/source degrees a={a} b={b}',
              zero(degree(other+[K,Y0])-expected(a,b,'unlike'))
              and zero(degree(other+[K,Y2])-expected(a,b,'equal')))
        check(f'fused L degree with explicit cutoff a={a} b={b}',
              zero(degree(other+[L],2*h)-expected(a,b,'unlike')))
    if b>=2:
        check(f'two-tail degree a={a} b={b}',
              zero(degree([A]*a+[B]*(b-2)+[K,K])-expected(a,b,'equal')))
    check(f'fused P degree a={a} b={b}',
          zero(degree(cur+[P],2*h)-expected(a,b,'unlike')))
    if a:
        check(f'even-root tail degree a={a} b={b}',
              zero(degree([A]*(a-1)+[B]*b+[FA])-expected(a,b,'even')))
check('source-only relative measure is two',zero(degree([Y0,Y2])-2*h)
      and degree([Y0,Y0])==4)

# Off-diagonal analytic family: independently compute the minority derivative.
z,w=S.symbols('z w')
aa=S.symbols('a1:4')
G=(z-w)**2/S.prod((w-a)**2 for a in aa)
logR=-rho*sum(1/(z-a) for a in aa)
f=sum(1/(z-a) for a in aa)-(2/rho)/(z-w)
check('even-root scalar derivative holds off g=rho',
      zero(G*f+(S.diff(G,z)+logR*G)/rho))
Si=[(z-a)**2/S.prod((a-b)**2 for b in aa if b!=a) for a in aa]
for a,si in zip(aa,Si):
    fi=sum(1/(z-b) for b in aa)-(2/rho)/(z-a)
    check(f'even-root tail derivative at {a}',zero(si*fi+(S.diff(si,z)+logR*si)/rho))
check('minority hole and moving-hole first variation vanish',
      zero(G.subs(w,z)) and zero(S.diff(G,z).subs(w,z)))

# Global source area finite part includes the infinity coefficient.
# For a rational holomorphic G with double poles on the real line,
# pi^-1 FP integral G = sum double-pole coefficients - [w^-2]_infinity.
aa_num=[S.Integer(-2),S.Integer(1),S.Integer(4)]
for count in [2,3]:
    gg=G.subs(dict(zip(aa,aa_num))) if count==3 else (z-w)**2/((w+2)**2*(w-1)**2)
    infinity=S.limit(w*w*gg,w,S.oo)
    check(f'even-root infinity term with {count} A currents',infinity==(1 if count==2 else 0))
gamma=t/(1-t)
check('combined even-root derivative coefficient',zero(1+gamma-1/(1-t)))
check('AA pair combined coefficient',zero(-2-2*gamma+2/(1-t)))
check('AA pair subtraction primitive',zero(S.diff(eps**(g-1)/(g-1),eps)-eps**(g-2)))
check('AAAB is universally marginal before derivative improvement',
      zero(expected(3,1,'even')-3*(g-rho)))
check('proper AAA positive on working patch',positive_patch(3*g-2))
check('proper AAB positive on working patch',positive_patch(g-2*rho+1))

delta=S.symbols('delta')
inner_primitive=R**(1+3*delta)/(1+3*delta)
check('nested integration gives improved inner denominator',
      zero(S.diff(inner_primitive,R)-R**(3*delta)))
outer=R**(3*delta+2*h-2)/((1+3*delta)*(3*delta+2*h-2))
check('nested integration has no residual delta pole',
      zero(S.diff(outer,R)-R**(2*h-4)*inner_primitive)
      and zero(S.limit(delta*outer,delta,0)))

# Check pole-times-zero danger directly, without trusting the shuffle output.
s,u=S.symbols('s u', nonzero=True)
c=1+s+1/s
# Labelled permutation sum, exchange weights from adjacent inversions.
from itertools import permutations
words=[('B','B','B','A'),('B','B','A','B'),('B','A','B','B'),('A','B','B','B')]
M=S.zeros(4)
for col,word in enumerate(words):
    for perm in permutations(range(4)):
        spatial=tuple(word[i] for i in perm)
        row=words.index(spatial)
        inv=[(perm[i],perm[j]) for i in range(4) for j in range(i+1,4) if perm[i]>perm[j]]
        factor=S.prod(u if word[i]==word[j] else 1/s for i,j in inv)
        M[row,col]+=factor
defect=M*S.Matrix([1,-c,c,-1])
Dstar=(s-1)*(s+1)*(s-u)*(s*s-u)*(u+1)/s**4
check('independent labelled shuffle defect',all(zero(x-y) for x,y in zip(defect,[0,-Dstar,Dstar,0])))
check('a chamber pole would survive the shuffle zero',
      zero(S.diff(Dstar,u).subs(u,s)*(-S.I*S.pi*s)/3
           -S.I*S.pi*(s-1/s)**2/3))
print(f'{checks} PASS, 0 FAIL')
