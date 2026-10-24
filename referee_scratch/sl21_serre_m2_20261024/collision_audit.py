"""Enumerate all current subsets and order-two collision power bounds.

Counts relative dimensions, not the dimension of a common current coordinate.
Checks a physical parameter patch by exact polynomial positivity, not sampling.
"""
import sympy as S
from itertools import permutations

passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    passed+=int(bool(ok));failed+=int(not bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)
def eq(name,a,b=0):check(name,S.expand(S.together(a-b).as_numer_denom()[0])==0)

b,t,g,rho,h=S.symbols('b t g rho h',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I]);a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I]);d=-(a0+a2)/b**2
dot=lambda x,y:S.simplify(x.dot(y)).subs(b*b,t)
eq('fresh root norm a0',dot(a0,a0),-1);eq('fresh root norm a2',dot(a2,a2),-1)
eq('fresh unlike screening dot product',dot(a0,a2),1-t)
eq('fresh dual-current norm',dot(d,d),-2/t)
eq('fresh dual/source pairing',dot(d,a0),1)

V={'A':(-1,0,0,0),'B':(1,0,0,1),'K':(1,0,1,0),'L':(1,1,1,1),
   'H':(-1,-1,-1,0),'0':(0,1,0,0),'2':(0,0,1,0),
   'P':(0,1,1,0),'V':(0,-1,-1,0)}
def degree(labels,sources,cutoff_power=0):
    vs=[V[x] for x in labels];deg=0
    for i,(qi,mi0,mi2,_) in enumerate(vs):
        for qj,mj0,mj2,_ in vs[i+1:]:
            screen=mi0*mj0+mi2*mj2+(h-1)*(mi0*mj2+mi2*mj0)
            base=rho*qi*qj
            if qi and qj and qi==qj:base+=g-rho
            deg+=base-qi*(mj0+mj2)-qj*(mi0+mi2)+2*screen
    currents=len(labels)-sources
    return S.expand(deg-sum(v[3] for v in vs)+(currents-1)+2*sources+cutoff_power)

subsets=sorted(set((a,b) for A,B in ((1,3),(3,1)) for a in range(A+1) for b in range(B+1) if a+b))
formulas={}
for a,b in subsets:
    base=g*(S.binomial(a,2)+S.binomial(b,2))-rho*a*b
    unlike=base+3*a-2*b+2*h+1
    equal=base+3*a-2*b+5
    even=base-a+2*b+1
    labels=['A']*a+['B']*b
    eq(f'unlike double-source count ({a},{b})',degree(labels+['0','2'],2),unlike)
    eq(f'equal double-source count ({a},{b})',degree(labels+['2','2'],2),equal)
    eq(f'even-root count ({a},{b})',degree(labels+['V'],1),even)
    eq(f'potential with cutoff count ({a},{b})',degree(labels+['P'],1,2*h),unlike)
    if b:
        replaced=['A']*a+['K']+['B']*(b-1)
        eq(f'tail/unlike-source count ({a},{b})',degree(replaced+['0'],1),unlike)
        eq(f'tail/equal-source count ({a},{b})',degree(replaced+['2'],1),equal)
        eq(f'current counterterm count ({a},{b})',degree(['A']*a+['L']+['B']*(b-1),0,2*h),unlike)
    if b>=2:eq(f'two-tail count ({a},{b})',degree(['A']*a+['K','K']+['B']*(b-2),0),equal)
    if a:eq(f'even-root tail count ({a},{b})',degree(['H']+['A']*(a-1)+['B']*b,0),even)
    formulas[(a,b)]=(unlike,equal,even)

# The zero-current source pair has two relative coordinates, NOT three.
eq('unlike source pair without a current',2*(h-1)+2,2*h)
eq('equal source pair without a current',2+2,4)

# Every new double-source degree is positive throughout 2<t<3. Multiplying
# by positive t reduces the question to a degree-two polynomial on [0,1].
x=S.symbols('x',real=True)
def bernstein_positive(expr):
    poly=S.Poly(S.cancel(t*expr.subs({g:2/t,rho:2/t,h:t})).subs(t,2+x),x)
    if poly.degree()>2:raise ValueError(poly)
    c0=poly.nth(0);c1=poly.nth(1);c2=poly.nth(2)
    B=[c0,c0+c1/2,c0+c1+c2]
    return all(v>=0 for v in B) and any(v>0 for v in B),B
for (a,b),(unlike,equal,even) in formulas.items():
    for name,value in [('unlike',unlike),('equal',equal)]:
        ok,B=bernstein_positive(value)
        check(f'positive new {name} face ({a},{b}) for 2<t<3',ok)
    if (a,b)==(1,0):eq('single A/even-root is its fixed Cauchy pole',even,0)
    elif (a,b)==(2,0):eq('new AA/even-root pair subtraction exponent',even,g-1)
    elif (a,b)==(3,1):eq('new marginal full AAAB/even-root face',even,3*(g-rho))
    else:
        ok,B=bernstein_positive(even)
        check(f'positive other even-root face ({a},{b}) for 2<t<3',ok)

# Check the uniform nested-face improvement once the inherited marginal
# first-order scalar has been written as a minority-coordinate derivative.
delta=S.symbols('delta')
inner_power=1+3*delta
outer_derivative_degree=2*h-5
outer_radial_measure=1
eq('nested first-order face has the full unlike joint degree',inner_power+outer_derivative_degree+outer_radial_measure+1,3*delta+2*h-2)
eq('equal-source nested face has positive joint degree',(inner_power+outer_derivative_degree+outer_radial_measure+1).subs(h,2),3*delta+2)

# Diagnostic of the change caused by an independent finite P perturbation.
eq('finite P would have a marginal SECOND jet in BBBA',degree(['B','B','B','A','P'],1)+2,3*(g-rho))

# Audit the permutation coefficient afresh from labelled inversions.
u,s=S.symbols('u s',nonzero=True)
words=['BBBA','BBAB','BABB','ABBB'];matrix=S.zeros(4)
for col,word in enumerate(words):
    for ordering in permutations(range(4)):
        spatial=''.join(word[i] for i in ordering)
        phase=1
        for left in range(4):
            for right in range(left+1,4):
                if ordering.index(left)>ordering.index(right):phase*=u if word[left]==word[right] else 1/s
        matrix[words.index(spatial),col]+=phase
co=S.Matrix([1,-1-s-1/s,1+s+1/s,-1])
dd=(s-1)*(s+1)*(s-u)*(s*s-u)*(u+1)/s**4
for i in range(4):eq(f'fresh exact two-exponent shuffle row {i}',(matrix*co)[i],[0,-dd,dd,0][i])
print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
