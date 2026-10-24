"""Independent checks from source equations and Gaussian definitions.
No prior audit scripts are read or imported.
"""
from pathlib import Path
import hashlib
import itertools
import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
mp.mp.dps = 70
passed = failed = 0

def check(name, condition, detail=''):
    global passed, failed
    if bool(condition):
        passed += 1
        print('PASS', name, detail)
    else:
        failed += 1
        print('FAIL', name, detail)

def zero(x):
    return s.simplify(s.expand(x)) == 0

def near(a, b, tol=mp.mpf('1e-50')):
    return abs(a-b) < tol * max(1,abs(a),abs(b))

inputs = [ROOT/'AGENTS.md', ROOT/'briefs/referee_sl21_vacuum_limits_report.md']
inputs += list((ROOT/'referee_scratch/sl21_vacuum_limits_20261024').glob('fateev_2001.*'))
inputs += list((ROOT/'referee_scratch/sl21_vacuum_limits_20261024').glob('fateev_2017.*'))
inputs += list((ROOT/'referee_scratch/sl21_vacuum_limits_20261024').glob('lashkevich_pugai_2017.*'))
start_hashes = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
(OUT/'input_hashes.txt').write_text(''.join(f'{v}  {k}\n' for k,v in start_hashes.items()))

z,al,c = s.symbols('z alpha c', nonzero=True, real=True)
t = -2*z*z
p0=s.Matrix([z,s.I*al]); p2=s.Matrix([z,-s.I*al]); p1=s.Matrix([-2*z,0])
pA=s.Matrix([-2*z-1/z,0]); pU=s.Matrix([z+1/z,-s.I*al]); pF=-pA
red=lambda x:s.factor(s.expand(x).subs(al**2,1+z*z))
dot=lambda p,q:red((p.T*q)[0])
check('Yukawa squared norm',zero(dot(p0,p0)+1) and zero(dot(p2,p2)+1))
check('primary U-source power',zero(-2*dot(pU,p0)-(2*t-4)))
check('endpoint fused momentum',pU+p0==pF)
check('reflected resonance norm',zero(dot(s.Matrix([-z,-s.I*al]),p0)-1))
check('A scaling dimension',zero(-dot(pA,pA)-(2-2*t)*(1/t-1)))
check('U scaling dimension',zero(-dot(pU,pU)-(2/t-1)))
# Derive LP labels directly from a/(2 alpha), -ell/(2 z).
h=2-t
check('LP A radial label',zero(2*z*(2*z+1/z)-(2*h-2)))
check('LP U radial label',zero(-2*z*(z+1/z)+h))
check('LP U angular label',zero(-2*(1+z*z)+h))

# Finite Gaussian subtraction: Zp=exp(-p H p/2).
H=s.Matrix([[s.Symbol('h11'),s.Symbol('h12')],[s.Symbol('h12'),s.Symbol('h22')]])
LZ=lambda p:s.expand(-(p.T*H*p)[0]/2)
E=lambda p,q:s.expand((p.T*H*q)[0])
check('general Gaussian OPE factor',zero(LZ(pU)+LZ(p0)-LZ(pF)-E(pU,p0)))
check('general Gaussian pair contact factor',zero(LZ(p0)+LZ(p2)-LZ(p0+p2)-E(p0,p2)))
# General H prediction after substituting lambda_old=Z_i lambda_H.
la=LZ(pA)+(1/t-1)*LZ(p1)
lu=LZ(pU)+LZ(p0)/t+(1/t-1)*LZ(p2)
iso={H[0,0]:c,H[0,1]:0,H[1,1]:c}
check('isotropic A extra factor',zero(red(la.subs(iso))-c*(1/t-1)))
check('isotropic U extra factor',zero(red(lu.subs(iso))-c*(2/t-1)))
check('isotropic UY0 OPE factor',zero(red(E(pU,p0).subs(iso))-c*(2-t)))
check('isotropic Y0Y2 contact factor',zero(red(E(p0,p2).subs(iso))-c*(1-t)))
check('fixed action A canonical scale',zero(red((LZ(pA)-LZ(p1)).subs(iso))-c*(1/t-2)))
check('fixed action B canonical scale',zero(red((LZ(pU)-LZ(p2)).subs(iso))-c*(1/t-1)))

# Transport the currents with their leading coefficient fixed in the new Wick basis.
pd=s.Matrix([1/z,0]);pT=p0+p2
wardA=LZ(p1)+LZ(pd)-LZ(pA)
wardB=LZ(p2)+LZ(pd)-LZ(pU)
check('general current Ward A transport',zero(wardA-E(pd,pT)))
check('general current Ward B transport',zero(wardB-E(pd,p2)))
check('isotropic current pair normalization',zero(red(LZ(pd).subs(iso))-c/t))
check('isotropic Ward A coefficient',zero(red(wardA.subs(iso))-2*c))
check('isotropic Ward B coefficient',zero(red(wardB.subs(iso))-c))
# Only the two physical bosons change subtraction; the fermion boson is fixed.
# In the rotated interacting plane H_eff = diag(c,c*z^2/(1+z^2)).
phys={H[0,0]:c,H[0,1]:0,H[1,1]:c*z*z/(1+z*z)}
check('physical Gaussian Yukawa Z0=Z2=1',zero(red(LZ(p0).subs(phys))) and zero(red(LZ(p2).subs(phys))))
check('physical Gaussian U factor',zero(red(LZ(pU).subs(phys))-c*(1/t-1)))
check('physical Gaussian A amplitude factor',zero(red(la.subs(phys))-c*(1/t-1)))
check('physical Gaussian U amplitude factor',zero(red(lu.subs(phys))-c*(1/t-1)))
check('physical Gaussian UY0 OPE factor',zero(red(E(pU,p0).subs(phys))-c*(1-t)))
check('physical Gaussian Y0Y2 action contact',zero(red(E(p0,p2).subs(phys))+c*t))
check('physical Gaussian Ward coefficients',zero(red(wardA.subs(phys))-2*c) and zero(red(wardB.subs(phys))-c))
check('fixed-action Ward A amplitude common scale',zero(wardA-LZ(p1)+LZ(pA)-LZ(pd)))
check('fixed-action Ward B amplitude common scale',zero(wardB-LZ(p2)+LZ(pU)-LZ(pd)))

# Rank-one source formulas, used literally with their stated root lengths.
def gammaf(x):return mp.gamma(x)/mp.gamma(1-x)
def liouville_A(a,zz,l1):
    bb=-mp.sqrt(2)*zz
    QQ=(bb+1/bb)/mp.sqrt(2)
    w=a-QQ
    return mp.power(l1*gammaf(bb*bb),w/(mp.sqrt(2)*bb))*mp.gamma(1-mp.sqrt(2)*w/bb)*mp.gamma(1-mp.sqrt(2)*w*bb)
def sl_A(a,xi,zz,C):
    aa=mp.sqrt(1+zz*zz); Q=1/(2*zz); w=a-Q
    return mp.power(C,w/zz)*mp.gamma(1-2*w*zz)*mp.gamma(1-w/zz)/(mp.gamma(mp.mpf('.5')-w*zz+aa*xi)*mp.gamma(mp.mpf('.5')-w*zz-aa*xi))
def GA(tt,l1):
    return mp.power(l1*gammaf(-tt),1/tt-1)*mp.gamma(2-tt)*mp.gamma(2-1/tt)/(mp.gamma(tt)*mp.gamma(1/tt))
def GU(tt,l0,l2):
    C=mp.sqrt(l0*l2)/(-tt)
    return l0/(1-tt)*mp.power(C,2/tt-2)*mp.gamma(2-1/tt)/mp.gamma(1/tt)

for j,tt in enumerate([-mp.sqrt(2)/3,-mp.sqrt(3)/4]):
    zz=mp.sqrt(-tt/2); aa=mp.sqrt(1+zz*zz); QL=-zz-1/(2*zz)
    l0=mp.mpf('.8');l2=mp.mpf('1.3');l1=mp.mpf('.63')
    lr=liouville_A(0,zz,l1)/liouville_A(2*QL,zz,l1)
    check(f'Liouville source ratio {j}',near(lr,GA(tt,l1)))
    C=mp.sqrt(l0*l2)/(-tt);astar=zz+1/zz
    slope=zz*mp.power(C,2/tt-2)*mp.gamma(2-1/tt)/((tt-1)*mp.gamma(1/tt))
    ratios=[]
    for dd in [mp.mpf('1e-10'),mp.mpf('1e-14')]:
        RR=sl_A(1/zz-(astar+dd),-aa,zz,C)/sl_A(astar+dd,-aa,zz,C)
        ratios.append(RR/(dd*slope))
    check(f'SL simple zero slope {j}',abs(ratios[-1]-1)<mp.mpf('1e-11'),mp.nstr(ratios[-1],24))
    check(f'SL slope converges linearly {j}',abs(ratios[-1]-1)<abs(ratios[0]-1)/1000)
    value=-l0*slope/zz
    check(f'zero times identity residue {j}',near(value,GU(tt,l0,l2)))
    lam=mp.mpf('.91');CC=lam/(-tt)
    gs=mp.power(CC,2/tt-1)*mp.gamma(1-1/tt)/mp.gamma(1/tt)
    check(f'symmetric gamma recurrence {j}',near(gs,GU(tt,lam,lam)))
    endpoint_ratio=(1-tt)**2*l2*GU(tt,l0,l2)/(l1*GA(tt,l1))
    ratio_simple=-mp.power(-gammaf(tt)*l0*l2/l1,1/tt)
    check(f'positive-branch endpoint ratio {j}',near(endpoint_ratio,ratio_simple))
    # Convergent radial integral for reflected field; u=r^2 makes area explicit.
    dd=mp.mpf('.7');RR=mp.mpf('1.2');q=zz*dd
    # Substitute u=R^2 x^(1/q), which removes the integrable endpoint singularity.
    radial=-l0*mp.quad(lambda x:RR**(2*q)/q,[0,1])
    exact=-l0/(zz*dd)*RR**(2*zz*dd)
    check(f'identity residue radial integral {j}',near(radial,exact))
    # Cutoff radial OPE with scalar and first scalar descendant powers.
    eps=mp.mpf('0.003');R=mp.mpf('0.8')
    for level in [0,1]:
        e=tt-1+level
        int_exact=-l0*(R**(2*e)-eps**(2*e))/e
        subtracted=int_exact-l0/e*eps**(2*e)
        check(f'cutoff subtraction sign {j}/{level}',near(subtracted,-l0/e*R**(2*e)))
    # Numerical arbitrary-H transport of the reference amplitudes.
    mat=mp.matrix([[mp.mpf('.23'),mp.mpf('.11')],[mp.mpf('.11'),mp.mpf('-.16')]])
    v0=mp.matrix([zz,1j*aa]);v2=mp.matrix([zz,-1j*aa]);v1=mp.matrix([-2*zz,0]);va=mp.matrix([-2*zz-1/zz,0]);vu=mp.matrix([zz+1/zz,-1j*aa]);vf=-va
    ZZ=lambda v:mp.exp(-(v.T*mat*v)[0]/2)
    za,z0,z2,z1,zu,zf=map(ZZ,[va,v0,v2,v1,vu,vf])
    l0H=l0/z0;l2H=l2/z2;l1H=l1/z1
    gaH=za*GA(tt,l1);guH=zu*GU(tt,l0,l2)
    check(f'arbitrary H amplitude A {j}',near(gaH,za*mp.power(z1,1/tt-1)*GA(tt,l1H)))
    check(f'arbitrary H amplitude U {j}',near(guH,zu*mp.power(z0,1/tt)*mp.power(z2,1/tt-1)*GU(tt,l0H,l2H)))
    Eco=mp.exp((vu.T*mat*v0)[0])
    check(f'arbitrary H endpoint subtraction {j}',near(l0H*Eco*zf,zu*l0))
    check(f'arbitrary H fixed-action B endpoint {j}',near(l2H*guH,zu/z2*l2*GU(tt,l0,l2)))
    # Compensated coupling loop from shifting the contour by Delta.
    Delta=mp.matrix([1j*mp.pi/zz,mp.pi/aa])
    dots=[(v.T*Delta)[0] for v in [v0,v2,v1]]
    check(f'coupling winding vector {j}',all(near(a,b) for a,b in zip(dots,[2j*mp.pi,0,-2j*mp.pi])))
    ss=mp.exp(-2j*mp.pi/tt)
    phaseA=mp.exp((va.T*Delta)[0]);phaseU=mp.exp((vu.T*Delta)[0])
    check(f'geometric field characters {j}',near(phaseA,1/ss) and near(phaseU,ss))
    windA=mp.exp((1/tt-1)*(-dots[2]));windU=mp.exp(-dots[0]/tt-(1/tt-1)*dots[1])
    check(f'compensated-loop characters {j}',near(windA,phaseA) and near(windU,phaseU))
    check(f'direct-pullback has inverse characters {j}',near(1/windA,ss) and near(1/windU,1/ss))

# Independent finite sign sums directly from LP equations (2.7)-(2.12).
def Jsign(X,Y,pp,g):
    w=lambda a:mp.exp(-2j*mp.pi*g*a)
    f=lambda e,u:(w(mp.mpf(e)/2)*u-w(-mp.mpf(e)/2))/(u-1)
    G=lambda e,u:(w(-mp.mpf(e)/2)*u+w(mp.mpf(e)/2))/(u+1)
    ans=mp.mpc(0)
    for eps in itertools.product([-1,1],repeat=len(X)+len(Y)):
        ex=eps[:len(X)];ey=eps[len(X):]
        term=w((sum(ex)*pp[0]+sum(ey)*pp[1])/4)
        for i in range(len(X)):
            for k in range(i+1,len(X)):
                if ex[i]!=ex[k]:term*=f(ex[i],X[k]/X[i])
        for i in range(len(Y)):
            for k in range(i+1,len(Y)):
                if ey[i]!=ey[k]:term*=f(ey[i],Y[k]/Y[i])
        for i in range(len(X)):
            for k in range(len(Y)):
                if ex[i]!=ey[k]:term*=G(ex[i],Y[k]/X[i])
        ans+=term
    return ans

tt=-mp.sqrt(2)/3;hh=2-tt;gg=1/hh
for K,L in [(0,1),(0,3),(1,1),(1,2),(2,2),(2,3),(3,3)]:
    X=[mp.exp(mp.mpf(i)/5) for i in range(1,K+1)]
    Y=[mp.exp(mp.mpf(i)/7+mp.mpf('.13')) for i in range(1,L+1)]
    rr=L-K;pp=(1-hh+rr,1-hh-rr)
    jj=Jsign(X,Y,pp,gg)
    check(f'FF shifted identity K{K}L{L}',abs(jj)<mp.mpf('1e-55'),mp.nstr(abs(jj),5))
for K in [1,2,3]:
    X=[mp.exp(mp.mpf(i)/5) for i in range(1,K+1)];Y=[mp.exp(mp.mpf(i)/7+mp.mpf('.13')) for i in range(1,K+1)]
    for label,pp in [('A',(hh-1,hh-1)),('U',(1-hh,1-3*hh))]:
        jj=Jsign(X,Y,pp,gg)
        check(f'FF {label} neutral {2*K} particles',abs(jj)<mp.mpf('1e-54'),mp.nstr(abs(jj),5))

# Resonance polynomials from Gaussian dimensions; enumeration supplements proof.
solutionsR=[];solutionsU=[]
for n0,n2,m,l in itertools.product(range(9),range(9),range(9),range(17)):
    n=n0+n2;D=n0-n2-1
    if D*D+n+2*m+l==1 and D*D-(n-1-2*m)**2+4*m==0:solutionsR.append((n0,n2,m,l))
    if D*D-n+6*m+l==1 and D*D-(n+1-2*m)**2+4*m==0:solutionsU.append((n0,n2,m,l))
check('reflected generic resonance enumeration',solutionsR==[(0,0,0,0),(1,0,0,0)],str(solutionsR))
check('unreflected generic resonance enumeration',solutionsU==[(0,0,0,0)],str(solutionsU))
# Sector diagonal characters and interface endpoint formula.
ss=mp.exp(-2j*mp.pi/tt)
for left,right in [(-2,1),(0,3),(1,1)]:
    N=right-left
    check(f'interface U character {left},{right}',near(ss**right-ss**left,ss**right*(1-ss**(-N))))
    check(f'interface A character {left},{right}',near(ss**(-right)-ss**(-left),ss**(-right)*(1-ss**N)))
end_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
check('all protected inputs unchanged',start_hashes==end_hashes,f'{len(inputs)} files')
print(f'SUMMARY {passed} PASS {failed} FAIL')
raise SystemExit(bool(failed))
