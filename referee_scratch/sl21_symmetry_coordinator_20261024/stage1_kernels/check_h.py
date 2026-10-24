"""Fresh independent h-regulator checks, no prior research imports."""
import sympy as S
import mpmath as M
import hashlib,json
from pathlib import Path
t,h,q,eps,R=S.symbols('t h q eps R',nonzero=True)
c,e,c1,e1=S.symbols('c e c1 e1')
eta=h-t;d=-(c+e)/t;d1=-(c1+e1)/t;beta=t-1
out=Path(__file__).resolve().parent;root=out.parents[2]
inputs=[root/'briefs/referee_sl21_mixed_regulator_report.md',root/'briefs/referee_sl21_current_renorm_report.md',root/'briefs/referee_sl21_symmetry_coordinator_plan_20261024.md',root/'briefs/referee_sl21_coord_stage0_audit_report.md']
inputs+=list((out.parent/'stage0_worker').iterdir())
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs if p.is_file()}
(out/'protected_inputs.json').write_text(json.dumps(hashes,indent=2)+'\n')
def ck(n,x):
    ok=bool(x);print(('PASS' if ok else 'FAIL')+' '+n)
    if not ok:raise AssertionError(n)
F0logprime=c+d+(h-2)/q
F2logprime=e+d+(h-2)/q
actualS0=c+d+(t-2)/q
actualS2=c+beta*d+(2-t)/q
ck('node0 actual residue minus derivative is −eta/q',S.simplify(actualS0-F0logprime+eta/q)==0)
ck('node2 actual residue plus derivative is eta/q',S.simplify(actualS2+F2logprime-eta/q)==0)
N=((c1-e1)+(c-e)*(c+e))/2
RP=c*(c+e)-t*((c1+e1)+(c+e)**2)/2
dL=c1+c*(c+e)*(1-1/t)
ck('fresh bosonized Ward polynomial N−RP/t=partial L',S.simplify(N-RP/t-dL)==0)
ck('regulator-dependent cutoff coefficient is t/h',S.simplify(1-eta/h-t/h)==0)
ck('full cutoff polynomial remainder factors eta',S.simplify(N-RP/h-t*dL/h-eta*N/h)==0)
ck('source-pair analytic radial primitive',S.simplify(S.diff(R**(2*h)/(2*h),R)-R**(2*h-1))==0)
ck('endpoint leading radial primitive',S.simplify(S.diff(R**(2*h-2)/(2*h-2),R)-R**(2*h-3))==0)
ck('h-improvement specializes to original',S.simplify((t/h).subs(h,t)-1)==0)
ck('h-P coefficient specializes to original',S.simplify((1/h).subs(h,t)-1/t)==0)
ck('h-endpoint coefficient specializes to original',S.simplify((1/(h-1)).subs(h,t)-1/(t-1))==0)
labels={'A':(-1,0,0,0),'B':(1,0,0,0),'D':(0,1,0,0),'0':(0,0,1,0),'2':(0,0,0,1),'K':(1,0,0,1),'Q':(0,1,0,1),'L':(1,0,1,1),'M':(0,1,1,1),'P':(0,0,1,1)}
def coc(n,m):return (-1)**((n[0]+n[1])*(m[2]+m[3]))*S.exp(-S.I*S.pi*(n[0]*m[1]-n[1]*m[0])/t)
def word(w):
    n=(0,0,0,0);res=1
    for key in w:
        m=labels[key];res*=coc(n,m);n=tuple(a+b for a,b in zip(n,m))
    return S.simplify(res)
for w,expected in [('ABD02',1),('ABQ0',-1),('AKD0',1),('ABM',1),('ALD',1),('ABDP',1)]:ck('ordered channel cocycle '+w,word(w)==expected)
# Pair exponents from a compact bilinear Gram matrix on {a0,a2}.
G=S.Matrix([[-1,1-t],[1-t,-1]])
a0=S.Matrix([1,0]);a2=S.Matrix([0,1]);T=a0+a2;dv=-T/t;kv=dv+T;zero=S.zeros(2,1)
dot=lambda x,y:S.simplify((x.T*G*y)[0])
mom={'A':(-dv,zero),'B':(dv,zero),'D':(zero,dv),'0':(a0,a0),'2':(a2,a2),'K':(dv+a2,a2),'Q':(a2,dv+a2),'L':(kv,T),'M':(T,kv),'P':(T,T)}
orders=['ABD02','ABQ0','AKD0','ABM','ALD','ABDP']
data={}
for w in orders:
    data[w]={}
    for i,x in enumerate(w):
        for y in w[i+1:]:
            powers=[S.simplify(-dot(mom[x][j],mom[y][j])) for j in [0,1]]
            if (x in '2KQ' and y=='0') or (x=='0' and y in '2KQ'):powers=[S.simplify(p+eta) for p in powers]
            if powers!=[0,0]:data[w][x+y]=[str(p) for p in powers]
(out/'channel_powers.json').write_text(json.dumps(data,indent=2)+'\n')
ck('double source shifted pair exponents',data['ABD02']['02']==['h - 1','h - 1'])
ck('bar-tail/source asymmetric shifted powers',data['ABQ0']['Q0']==['h - 1','h - 2'])
ck('tail/source asymmetric shifted powers',data['AKD0']['K0']==['h - 2','h - 1'])
ck('A-fused endpoint exponent',dot(dv,kv)==2-2/t)
# A local Gaussian test of the fully subtracted eta anomaly coefficient.
# The second Taylor coefficient is N=1; replace its smooth radial test by exp(-r²).
# Its continued integral is Gamma(h), whereas raw integral converges only h>0.
M.mp.dps=50
def radial_fp(hv):
    # subtract constant near zero for −1<Re h<0; finite radial convention at R=1
    return 2*M.quad(lambda r:r**(2*hv-1)*M.expm1(-r*r),[0,1])+2*M.quad(lambda r:r**(2*hv-1)*M.exp(-r*r),[1,M.inf])+1/hv
for tv in [-M.sqrt(2)/2,-1/M.sqrt(3)]:
    ck('independent Mellin finite part equals Gamma at '+str(float(tv)),abs(radial_fp(tv)-M.gamma(tv))<M.mpf('1e-28'))
    for displacement in [M.mpf('.01'),M.mpf('.001')]:
        hv=tv+displacement;val=displacement*radial_fp(hv)
        print('INFO eta anomaly',float(tv),str(displacement),M.nstr(val,18))
        ck('anomaly divided by eta tends to regular Gamma',abs(val/displacement-M.gamma(tv))<50*displacement)
# Pole divisors from collision sectors: coefficient of Taylor order n.
n=S.symbols('n',integer=True,nonnegative=True)
degrees={'02':2*h,'B02':2*h-1,'BD02':2*h-3,'A02':2*h+4,'AB02':2*h+2-2/t,'AD02':2*h+2,'ABD02':2*h-2/t}
for key,degree in degrees.items():
    restricted=S.expand((degree+n).subs(h,t))
    ck('no Mellin divisor identically contains diagonal: '+key,S.diff(restricted,t)!=0)
print('INFO Pole enumeration establishes generic diagonal restriction, not absence of isolated resonance values.')
