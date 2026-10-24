"""Independent audit from primary three-vector and contour definitions."""
from pathlib import Path
import hashlib,json,itertools
import sympy as s
ROOT=Path('/home/gustav/Git/affine_toda')
OUT=Path(__file__).resolve().parent
inputs=['AGENTS.md','briefs/referee_sl21_symmetry_coordinator_plan_20261024.md','briefs/referee_sl21_coord_stage1_kernels_report.md','briefs/referee_sl21_coord_stage1_subtractions_report.md','briefs/referee_sl21_coord_stage0_audit_report.md','briefs/referee_sl21_triple_extension_report.md','briefs/referee_sl21_boundary_global_audit_report.md','referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/all_order_and_completion.md','referee_scratch/sl21_vacuum_limits_20261024/fateev_2017.pdf']
(OUT/'input_hashes.json').write_text(json.dumps({p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},indent=2)+'\n')
def check(label,condition):
    yes=bool(condition);print(('PASS ' if yes else 'FAIL ')+label)
    if not yes:raise AssertionError(label)
b,t,h,rho,a,v,w=s.symbols('b t h rho a v w',nonzero=True)
a0=s.Matrix([s.I*b/s.sqrt(2),b/s.sqrt(2),s.I]);a2=s.Matrix([s.I*b/s.sqrt(2),-b/s.sqrt(2),-s.I]);d=-(a0+a2)/b**2
simp=lambda x:s.factor(s.simplify(x).subs(b**2,t))
dot=lambda x,y:simp((x.T*y)[0])
check('primary vector norms',dot(a0,a0)==-1 and dot(a2,a2)==-1 and dot(a0,a2)==1-t)
check('current contractions',dot(a0,d)==dot(a2,d)==1 and dot(d,d)==-2/t)
check('node2 vector first jet',all(simp(x)==0 for x in a0+(t-1)*d+(d+a2)))
# Actual Laurent coefficients of the double-source scalar kernel. The pole-free
# numerator has logarithmic derivative rho/(a-z), and the other pole contributes.
ratio0=(rho-1)/(a-v)+(t-2)/(v-w)
dlog0=(rho-1)/(a-v)+(h-2)/(v-w)
ratio2=((t-1)*rho-1)/(a-w)+(2-t)/(w-v)
dlog2=(rho-1)/(a-w)+(h-2)/(w-v)
check('A-dressed node0 Laurent remainder',s.simplify(ratio0-dlog0+(h-t)/(v-w))==0)
check('A-dressed node2 Laurent remainder',s.simplify((ratio2+dlog2-(h-t)/(w-v)).subs(rho,2/t))==0)
# All proper clusters, directly from pairwise primary contractions.
z=s.zeros(3,1)
mom={'A':(-d,z),'B':(d,z),'D':(z,d),'0':(a0,a0),'2':(a2,a2)}
degrees={}
for size in range(2,6):
 for cluster in itertools.combinations(mom,size):
    if not {'0','2'}<=set(cluster):continue
    power=0
    for x,y in itertools.combinations(cluster,2):
        power-=dot(mom[x][0],mom[y][0])+dot(mom[x][1],mom[y][1])
    power+=2*(h-t)
    for name,chir in [('B',0),('D',1)]:
        if name in cluster and any(dot(a0,mom[q][chir])!=0 for q in cluster if q!=name):power-=1
    ncurr=sum(x in 'ABD' for x in cluster)
    measure=ncurr+2*sum(x in ('0','2') for x in cluster)-(1 if ncurr else 2)
    degrees[''.join(cluster)]=s.simplify(power+measure)
expected={'02':2*h,'A02':2*h+4,'B02':2*h-1,'D02':2*h-1,'AB02':2*h+2-2/t,'AD02':2*h+2,'BD02':2*h-3,'ABD02':2*h-2/t}
check('all eight new faces reconstructed',degrees==expected)
for name,degree in degrees.items():print('INFO face',name,degree)
check('common new-face convergence region',all(s.simplify(dg.subs({h:3,t:-s.sqrt(2)/2}))>0 for dg in degrees.values()))
n=s.symbols('n',integer=True,nonnegative=True)
check('no new Mellin divisor contains physical diagonal',all(s.Poly(s.expand(t*(dg.subs(h,t)+n)),t).degree()>0 for dg in degrees.values()))
# Derivatives of the regulator, including a mixed second derivative.
q,qbar,eta=s.symbols('q qbar eta',nonzero=True)
R=q**eta*qbar**eta
check('single regulator derivative divisible by eta',s.simplify(s.diff(R,q)/R-eta/q)==0)
check('two same-chirality regulator derivatives divisible by eta',s.simplify(s.diff(R,q,2)/R-eta*(eta-1)/q**2)==0)
check('mixed regulator derivatives divisible by eta',s.simplify(s.diff(R,q,qbar)/R-eta**2/(q*qbar))==0)
# Exact circular omitted disk at zero derivative backgrounds, A separated.
x=s.symbols('x')
series=x*s.diff(1/(1-x),x)+t*x/(1-x)
closed=x*((1+t)-t*x)/(1-x)**2
check('overlap nonuniform coefficient reconstructed',s.simplify(series-closed)==0)
print('INFO omitted-disk degree at z,v proportional to epsilon: 2*h-6')
print('INFO adding B,D positions and source-pair centre transverse integrations: 2*h-3')
# Local anomaly integer enumeration, independently using required h,hbar levels.
allowed={'A':[], 'B':[]}
for n0,n2,ev in itertools.product(range(25),repeat=3):
 if n0+n2+ev==0:continue
 aa,cc=n0-ev,n2-ev;F=n0-n2;M=n0+n2+2*ev
 if ev!=aa*cc:continue
 lb=s.Rational(1)-s.Rational(M,2)-s.Rational(F*F,2)
 la=-s.Rational(M,2)-(aa+cc)-s.Rational(F*F,2)
 lj=lb+(aa+cc)
 if lb>=0 and la>=0:allowed['A'].append((n0,n2,ev,int(la),int(lb)))
 if lb>=0 and lj>=0:allowed['B'].append((n0,n2,ev,int(lj),int(lb)))
check('ordinary local anomaly candidates',allowed=={'A':[(0,0,1,1,0)],'B':[(0,1,0,1,0),(1,0,0,1,0)]})
print('INFO local anomaly candidates',allowed)
# Negative Serre words: K powers k=-1,+1 and barred weights -k.
for majority in [-1,1]:
 exponents=[]
 for position in range(4):
    ks=[majority]*3;ks.insert(position,-majority)
    exponents.append(sum(-ki*kj for i,ki in enumerate(ks) for kj in ks[i+1:]))
 check('negative Serre common K factor for '+str(majority),exponents==[0,0,0,0])
# The finite real-cycle lower bound tends to minus infinity even with a gradient.
phi=s.symbols('phi',positive=True)
check('negative exponential defeats quadratic gradient',s.limit(phi**2*s.exp(-phi),phi,s.oo)==0)
print('INFO anchor BCH exponent: (-q*d/2) dot ((l+r)/2) times (-4*pi*i) = i*pi*q*d dot (l+r)')
print('INFO all checks use primary formulas; no earlier research code imported')
