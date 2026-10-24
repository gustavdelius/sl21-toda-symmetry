"""Fresh same-chirality Wick and permutation checks; no prior code imported."""
from collections import defaultdict
from itertools import permutations, combinations
import math
import cmath
import numpy as np
import sympy as S
from scipy.integrate import quad
from scipy.special import beta as beta_function

npas = nfail = 0
def check(name, value, detail=''):
    global npas,nfail
    good=bool(value)
    npas+=int(good); nfail+=int(not good)
    print(('PASS ' if good else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)
def eq(name,a,b=0):
    # Clear factored denominators before expansion. A multivariate gcd here
    # is unnecessary for an identity test and was prohibitively slow.
    numerator=S.together(a-b).as_numer_denom()[0]
    check(name,S.expand(numerator)==0)

t,rho,s=S.symbols('t rho s', nonzero=True)
b=S.symbols('b',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I])
a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I])
d=-(a0+a2)/b**2
def dot(a,bv):return S.simplify(a.dot(bv)).subs(b**2,t)
eq('a0 square',dot(a0,a0),-1)
eq('a2 square',dot(a2,a2),-1)
eq('a0 a2',dot(a0,a2),1-t)
eq('a0 d',dot(a0,d),1)
eq('d square',dot(d,d),-2/t)

# Pure current ordered integrals: sum all labelled permutations, with the
# canonical separated braid for every inversion. Identical species are kept
# labelled until their contributions are collected.
def symmetrize(word):
    out=defaultdict(lambda:S.Integer(0))
    for p in permutations(range(len(word))):
        power=0
        for i in range(len(p)):
            for j in range(i+1,len(p)):
                if p[i]>p[j]:
                    power+=1 if word[p[i]]==word[p[j]] else -1
        out[''.join(word[i] for i in p)]+=s**power
    return {w:S.expand(c) for w,c in out.items()}

matrices={}
for n in range(2,5):
    for na in range(n+1):
        words=sorted(set(''.join(p) for p in permutations('A'*na+'B'*(n-na))))
        M=S.Matrix([[symmetrize(w).get(v,0) for w in words] for v in words])
        matrices[n,na]=(words,M)
        # Exact generic specialization establishes full rank unless an
        # identically verified symbolic kernel is present.
        target=len(words)-(1 if n==4 and na in (1,3) else 0)
        check(f'permutation matrix generic rank n={n}, nA={na}',M.subs(s,2).rank()==target,
              f'rank={target}/{len(words)}')
        if n==4 and na in (1,3):
            repeat='A' if na==3 else 'B'; other='B' if na==3 else 'A'
            order=[repeat*3+other,repeat*2+other+repeat,repeat+other+repeat*2,other+repeat*3]
            coeffs=[1,-(s+1+1/s),s+1+1/s,-1]
            vec=S.Matrix([coeffs[order.index(w)] for w in words])
            check('quartic kernel '+repeat,all(S.cancel(x)==0 for x in M*vec))
            derived=M.nullspace()[0]
            derived=derived/derived[words.index(order[0])]
            check('quartic coefficients extracted from nullspace '+repeat,
                  all(S.cancel(derived[i]-vec[i])==0 for i in range(len(words))))
            print('DATA quartic '+repeat+' '+str(dict(zip(order,coeffs))),flush=True)

check('pure B is even, no generic quadratic truncation',symmetrize('BB')['BB'].subs(s,2)==3)
rn=-math.sqrt(2)/4
constant_jet_pair=quad(lambda r:(1-r)*r**rn,0,1,epsabs=1e-10)[0]
check('nonzero constant-jet B-square ordered integral',
      abs(constant_jet_pair-1/((rn+1)*(rn+2)))<1e-9 and constant_jet_pair>0)

# Derivative Wick products from monomer/dimer matchings. The derivative field
# has self-contraction +1/(z_i-z_j)^2; its contractions with B exponents give
# -1/(z_i-z_j), and with A exponents give +1/(z_i-a).
def wick_prefactor(z,external,C=lambda z:0):
    shifts=[C(x)+sum(q/(x-y) for y,q in external)
            -sum(1/(x-y) for j,y in enumerate(z) if i!=j) for i,x in enumerate(z)]
    def rec(ids):
        if not ids:return S.Integer(1)
        i,*rest=ids
        ans=shifts[i]*rec(rest)
        for j in rest:
            ans+=rec([k for k in rest if k!=j])/(z[i]-z[j])**2
        return ans
    return rec(list(range(len(z))))

z=S.symbols('z1:4'); a,w=S.symbols('a w'); c,cp=S.symbols('c cp')
F=S.symbols('F1:4')
two=wick_prefactor(z[:2],[],lambda x:F[z.index(x)])
eq('exact BB double pole cancels',two,F[0]*F[1]+(F[0]-F[1])/(z[0]-z[1]))
three=wick_prefactor(z,[],lambda x:F[z.index(x)])
dd2=sum(F[i]/S.prod(z[i]-z[j] for j in range(3) if j!=i) for i in range(3))
bounded=F[0]*F[1]*F[2]+sum(F[k]*(F[i]-F[j])/(z[i]-z[j]) for i,j,k in [(0,1,2),(0,2,1),(1,2,0)])+2*dd2
eq('BBB prefactor is a bounded divided-difference expression',three,bounded)
eq('BBBA with zero regular jets has no fully contracted term',wick_prefactor(z,[(a,1)]))
eq('BBBA constant regular jet starts at degree -1',wick_prefactor(z,[(a,1)],lambda x:c),
   c**3+c**2*sum(1/(x-a) for x in z))
eq('BBA constant regular jet',wick_prefactor(z[:2],[(a,1)],lambda x:c),
   c**2+c*sum(1/(x-a) for x in z[:2]))

# The actual UV convergence inequalities, after the cancellations just checked.
rg=-math.sqrt(2)/4
degrees={'AA/BB pair':rho+1,'AB pair':-rho,'AAA/BBB':3*rho+2,
         'AAB/BBA':1-rho,'AAAB/BBBA':S.Integer(2)}
for name,p in degrees.items():
    check('UV convergence '+name,float(p.subs(rho,rg))>0,str(p))

# Same-chirality B-tail leading singularity is not softened like BB.
eq('B K2 vertex power',-dot(d,d+a2),2/t-1)
eq('B K2 derivative contraction',-dot(a0,d+a2),t-2)
real_rho=S.symbols('real_rho',real=True)
check('no common absolute domain for AB and B-tail',S.reduce_inequalities([real_rho<0,real_rho>1],real_rho)==S.false,
      'AB needs Re(rho)<0, B-tail needs Re(rho)>1')

beta=S.symbols('beta')
for n in (2,3):
    expr=wick_prefactor(z[:n],[(w,beta)])
    expected=S.prod(beta-j for j in range(n))/S.prod(x-w for x in z[:n])
    eq(f'{n} B currents with one Yukawa, full contraction',expr,expected)

# Independent area integration of the two-double-pole kernel. Green's theorem
# predicts 2*pi/r^2 - 6*pi*epsilon^2/r^4 for two circular holes.
u,v=.5,-.5
rr=u-v
outer=3.
def primitive(x,y):
    zz=complex(x,y)
    return -(1/(zz-u)+1/(zz-v))/rr**2-2*(cmath.log(zz-u)-cmath.log(zz-v))/rr**3
def area_integral(e):
    def horizontal(y):
        xmax=math.sqrt(max(0.,outer*outer-y*y))
        if y<e:
            half=math.sqrt(e*e-y*y)
            intervals=[(-xmax,v-half),(v+half,u-half),(u+half,xmax)]
        else:intervals=[(-xmax,xmax)]
        return sum((primitive(hi,y)-primitive(lo,y)).real for lo,hi in intervals)
    return 2*(quad(horizontal,0,e,epsabs=2e-9,epsrel=2e-10,limit=200)[0]
              +quad(horizontal,e,outer,epsabs=2e-9,epsrel=2e-10,limit=200)[0])
for e in (.1,.05,.025):
    num=area_integral(e); exact=2*math.pi/rr**2-6*math.pi*e*e/rr**4
    check(f'direct source area integral, hole={e}',abs(num-exact)<2e-7,
          f'numeric={num:.12g}, boundary formula={exact:.12g}')

# Assemble signs from the given charge-cover cocycle and all current tails.
for n in (2,3):
    fall=S.prod(t-2-j for j in range(n-1))
    action=(-1)**(n+1)*(t-1)*fall
    tails=(-1)**n*t*fall
    eq(f'combined n={n} B product leading first-order coefficient',action+tails,(-1)**n*fall)
eq('BB first-order leading coefficient per r^(rho-2)',2*(t*(t-2)-(t-1)*(t-2)),2*(t-2))
eq('ordered BB contact subtraction simplifies to -2 t',
   (2*(t-2)/(rho-1)).subs(rho,2/t),-2*t)
rg=-math.sqrt(2)/4; tg=2/rg
R=.7
target=-2*tg*R**(rg-1)
for cutoff in (.07,.035):
    raw=2*(tg-2)*quad(lambda x:x**(rg-2),cutoff,R,epsabs=1e-8)[0]
    renormalized=raw-2*tg*cutoff**(rg-1)
    check(f'ordered BB leading subtraction epsilon={cutoff}',abs(renormalized-target)<1e-8)

u_gap,v_gap=S.symbols('u_gap v_gap',positive=True)
Dsum=1/(u_gap**2*(u_gap+v_gap)**2)+1/(u_gap**2*v_gap**2)+1/(v_gap**2*(u_gap+v_gap)**2)
eq('BBB source-first angular kernel',Dsum,2*(u_gap**2+u_gap*v_gap+v_gap**2)/(u_gap*v_gap*(u_gap+v_gap))**2)
rp=math.sqrt(2)
ang=quad(lambda x:(x*(1-x))**(rp-2)*(1-x+x*x),0,1,epsabs=2e-9,limit=200)[0]
ang_exact=(3*rp-1)/(2*(2*rp-1))*beta_function(rp-1,rp-1)
check('BBB angular integral independent quadrature',abs(ang-ang_exact)<2e-8,
      f'numerical={ang:.12g}, beta-function={ang_exact:.12g}')

# The primitive four-current/one-Yukawa coefficient is nonzero for node 2,
# but has a special derivative in the A coordinate when rho*t=2.
# Keeping all points symbolic provides a stronger check than a scale count.
p3=wick_prefactor(z,[(a,1),(w,beta)])
compact=beta*(beta-1)/S.prod(x-w for x in z)*(beta+1+(a-w)*sum(1/(x-a) for x in z))
eq('three B plus A plus Yukawa prefactor',p3,compact)
eq('node 0 primitive scalar vanishes',p3.subs(beta,1))
sample=S.factor(p3.subs({z[0]:1,z[1]:2,z[2]:3,a:0,w:4}))
eq('node 2 primitive scalar is not identically zero',sample,-beta*(beta-1)*(3*beta-19)/18)
print('DATA critical source prefactor example',sample,flush=True)

# Divide the derivative identity by its common exponential factor; do not
# differentiate fractional powers with SymPy's implicit branch assumptions.
sumv=sum(1/(x-a) for x in z)
kernel_without_R=(a-w)*compact.subs(beta,t-1)/S.prod(x-w for x in z)
derivative_without_R=(t-1)*(t-2)/rho*(2*(a-w)+rho*(a-w)**2*sumv)/S.prod((x-w)**2 for x in z)
eq('critical action kernel is an A derivative at rho*t=2',
   (kernel_without_R-derivative_without_R).subs(rho,2/t))

for i in range(3):
    rest=[z[j] for j in range(3) if i!=j]
    tail_pref=wick_prefactor(rest,[(a,1),(z[i],t-2)])
    D=1/S.prod((z[i]-x)**2 for x in rest)
    # Includes the local cocycle and the -t tail coefficient, as derived by
    # counting preceding/following current factors.
    tail=t*(a-z[i])*tail_pref/S.prod(x-z[i] for x in rest)
    deriv=t*(t-2)/rho*D*(2*(a-z[i])+rho*(a-z[i])**2*sumv)
    eq(f'critical tail {i} is an A derivative at rho*t=2',(tail-deriv).subs(rho,2/t))
eq('combined critical real-axis coefficient',
   (t*(t-2)-(t-1)*(t-2))/rho,(t-2)/rho)

# Numerical test of the source-first critical kernel with a real A point,
# using residues at the three B positions and direct differentiation.
def rational_residues(zz,aa):
    doubles=[]; residues=[]
    for i,zi in enumerate(zz):
        D=1/np.prod([(zi-zj)**2 for j,zj in enumerate(zz) if j!=i])
        doubles.append((aa-zi)**2*D)
        residues.append(-2*(aa-zi)*D-2*(aa-zi)**2*D*sum(1/(zi-zj) for j,zj in enumerate(zz) if j!=i))
    return np.array(doubles),np.array(residues)
for zz in ([.2,.6,1.1],[.2+.03j,.6+.02j,1.1+.01j]):
    zz=np.array(zz,dtype=complex); aa=-.7+.04j
    D,Res=rational_residues(zz,aa)
    check('source residue identity with '+('complex' if np.max(abs(zz.imag)) else 'real')+' current positions',
          abs(np.dot(zz,Res)+sum(D))<2e-10)
    exact=-math.pi*np.dot(np.conj(zz),Res)
    decomposition=math.pi*sum(D)+2j*math.pi*np.dot(zz.imag,Res)
    check('explicit contour-height term in source integral',abs(exact-decomposition)<2e-10)

print(f'RESULT {npas} PASS, {nfail} FAIL',flush=True)
raise SystemExit(bool(nfail))
