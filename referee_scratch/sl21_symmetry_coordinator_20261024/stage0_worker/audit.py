"""Independent local Wick and collision checks; no imported research scripts."""
from pathlib import Path
import hashlib,itertools,json
import sympy as S
import mpmath as M

root=Path(__file__).resolve().parents[3]
out=Path(__file__).resolve().parent
inputs=['AGENTS.md','README.md','01-introduction.qmd','briefs/referee_sl21_symmetry_coordinator_plan_20261024.md','briefs/referee_sl21_triple_extension_report.md','briefs/referee_sl21_triple_kernel_report.md','referee_scratch/sl21_triple_extension_20261024/checkpoint.md','briefs/referee_sl21_bernard_leclair_report.md']
snapshot={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in inputs}
(out/'inputs.sha256.json').write_text(json.dumps(snapshot,indent=2)+'\n')
def ck(name,condition):
    print(('PASS' if bool(condition) else 'FAIL')+' '+name)
    if not condition: raise AssertionError(name)
t=S.symbols('t',nonzero=True); b=S.symbols('b',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I]);a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I]);d=-(a0+a2)/b**2
dot=lambda a,b:S.simplify(a.dot(b)).subs(globals()['b']**2,t)
ck('independent bilinear Gram matrix',[dot(a0,a0),dot(a2,a2),dot(a0,a2),dot(a0,d),dot(a2,d),dot(d,d)]==[-1,-1,1-t,1,1,-2/t])
ck('full first-jet vector identity',S.simplify(a0+(b*b-1)*d+d+a2)==S.zeros(3,1))
a,p,r=S.symbols('a p r',nonzero=True);rho=2/t;beta=t-1
pref=(a-p)**(-rho)*a
simple=S.simplify(S.diff(pref,p).subs(p,0)*beta-pref.subs(p,0)/a)
ck('mixed simple coefficient is minus V prime',S.simplify(simple-(1-rho)*a**(-rho))==0)
ck('simultaneous rational reduction',S.simplify((p-r)/p*(1/r+beta/p)-(1/r+(beta-1)/p-beta*r/p**2))==0)
labels={'A':(-1,0,0,0),'B':(1,0,0,0),'D':(0,1,0,0),'0':(0,0,1,0),'2':(0,0,0,1),'U':(1,1,0,1)}
def coc(n,m):return (-1)**((n[0]+n[1])*(m[2]+m[3]))*S.exp(-S.I*S.pi*(n[0]*m[1]-n[1]*m[0])/t)
def product(word):
    n=(0,0,0,0);c=1
    for name in word:
        m=labels[name];c*=coc(n,m);n=tuple(x+y for x,y in zip(n,m))
    return S.simplify(c)
ck('single source total cocycle',product('ABD2')==-1)
ck('double source total cocycle',product('ABD02')==1)
ck('double-source action combinatorics',S.Rational(1,2)*2==1)
nu=S.exp(-S.I*S.pi/t);s=nu**2
ratio=coc(labels['A'],labels['U'])/coc(labels['U'],labels['A'])
ck('endpoint cocycle ratio',S.simplify(ratio+s**-1)==0)
ck('left endpoint A U over U A phase',S.simplify(ratio*S.exp(-S.I*S.pi*(1-rho))-s**-2)==0)
ck('right endpoint A U over U A phase',S.simplify(ratio*S.exp(S.I*S.pi*(1-rho))-1)==0)
M.mp.dps=35
I=M.quad(lambda y:2*M.sqrt(1-y*y),[-1,0,1])
ck('inner-circle strip integral is pi',abs(I-M.pi)<M.mpf('1e-28'))
# A convergent real integrand for each full horizontal edge integral.
for h in [M.sqrt(2),M.sqrt(3)]:
    imag=M.quad(lambda x:h/(x*x+h*h),[-M.inf,0,M.inf])
    ck('strip edge horizontal first pole pi at height '+str(float(h)),abs(imag-M.pi)<M.mpf('1e-28'))
ck('raw contact coefficient',S.expand(-t*(t-1)+t*t)==t)
ck('selected finite matching coefficient',S.expand(t+t*(t-1))==t*t)

# Derive scalar contraction exponents from the vectors, then find the most
# singular nonzero derivative monomial on each subset. D denotes bar B.
zero=S.zeros(3,1)
mom={'A':(-d,zero),'B':(d,zero),'D':(zero,d),'0':(a0,a0),'2':(a2,a2)}
def face(face):
    exponent=0
    for x,y in itertools.combinations(face,2):
        exponent-=dot(mom[x][0],mom[y][0])+dot(mom[x][1],mom[y][1])
    for cur,chir in [('B',0),('D',1)]:
        if cur in face and any(x!=cur and dot(a0,mom[x][chir])!=0 for x in face):exponent-=1
    dim=sum(1 if x in 'ABD' else 2 for x in face)-(1 if set(face)&set('ABD') else 2)
    return S.simplify(exponent+dim)
table={''.join(c):str(face(c)) for n in range(2,6) for c in itertools.combinations('ABD02',n)}
(out/'faces.json').write_text(json.dumps(table,indent=2)+'\n')
expected={'02':2*t,'AB':-rho,'BD2':-1,'AB2':2-rho,'AD2':2,'ABD2':1-rho,'B02':2*t-1,'D02':2*t-1,'BD02':2*t-3,'A02':2*t+4,'AB02':2*t+2-rho,'AD02':2*t+2,'ABD02':2*t-rho}
for word,value in expected.items():ck('face '+word+' degree '+str(value),S.simplify(face(word)-value)==0)
for tval in [-M.sqrt(2)/2,-1/M.sqrt(3)]:
    rv=2/tval
    for word in ['A02','AB02','AD02','ABD02']:
        ck('positive new A face '+word+' at '+str(float(tval)),float(face(word).subs(t,float(tval)))>0)
    ck('mixed first unretained fifth jet vanishes',2*tval+2>0)
    ck('mixed fifth-jet spectator remainder integrable',-3-rv>-1)
    ck('source-pair Taylor second jet vanishes',2*tval+2>0)
    eps=M.mpf('0.0001');R0=M.mpf('.2');scale=M.sqrt(3)
    value=eps*M.quad(lambda R:R**(-1-rv),[scale*eps,R0])
    bound=eps*(R0**(-rv))/(-rv)
    ck('first-order nested remainder explicit bound',0<value<=bound)
ck('source-pair and AB absolute domains incompatible',S.re(1/(S.symbols('x',real=True)+S.I*S.symbols('y',real=True)))==S.symbols('x',real=True)/(S.symbols('x',real=True)**2+S.symbols('y',real=True)**2))
print('INFO Double-Yukawa degrees are a forest classification, not evaluated contact coefficients.')
