"""Independent resonance and endpoint matching checks from primary formulas.

These calculations do not prove the massive reflection identity or identify
the original Toda Hilbert space with a sine-Liouville/CSG completion.
"""
import sympy as S
import mpmath as mp
mp.mp.dps = 90
count = 0

def check(name, ok):
    global count
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        raise AssertionError(name)
    count += 1

def zero(name, expr):
    check(name, S.simplify(expr) == 0)

t,z,e,r,R,lam = S.symbols('t z e r R lam', nonzero=True)
alpha2=1+z*z
astar=z+1/z
Q=1/(2*z)
partner=2*Q-(astar+e)
zero('SL reflection partner', partner-(-z-e))
zero('opposite screening contraction', partner*z+alpha2-(1-z*e))
zero('resonant radial antiderivative', S.diff(r**(2*z*e)/(2*z*e),r)-r**(-1+2*z*e))
zero('resonant source simple-pole residue', S.limit(e*(-lam*R**(2*z*e)/(z*e)),e,0)+lam/z)

# The source pair is normal ordered at r=0: only the angular scalar terms
# 1 and r*bar(r) C*bar(C) contribute at the displayed orders.
zero('U-Y0 contraction exponent', S.expand(astar*z+alpha2).subs(z*z,-t/2)-(2-t))
zero('leading angular radial primitive', S.diff(r**(2*t-2)/(t-1),r)-2*r**(2*t-3))
zero('descendant angular radial primitive', S.diff(r**(2*t)/t,r)-2*r**(2*t-1))
u,F,CBF,l0,l2,cut = S.symbols('U F CbarCF lambda0 lambda2 cutoff', nonzero=True)
Umin = u-l0*cut**(2*t-2)*F/(t-1)-l0*cut**(2*t)*CBF/t
BB = t*t*l2*u-t*t*l0*l2*cut**(2*t-2)*F/(t-1)-t*l0*l2*cut**(2*t)*CBF
zero('canonical BB is t squared lambda2 times minimally subtracted U', BB-t*t*l2*Umin)

# Generic irrational coupling: independent integer constant/v coefficients.
v=S.symbols('v')
n0,n2,m,level=S.symbols('n0 n2 m level', integer=True, nonnegative=True)
n=n0+n2
D=n0-n2-1
rad=n-1-2*m
Delta_fused=(1+v)*D*D-v*rad*rad
Delta_couplings=n+m*(2+4*v)
resonance=S.expand(Delta_fused+Delta_couplings+level-1)
zero('reflected resonance constant coefficient', resonance.coeff(v,0)-(D*D+n+2*m+level-1))
solutions=[]
for i in range(7):
    for j in range(7):
        for k in range(4):
            for ell in range(4):
                d=i-j-1; a=i+j-1-2*k
                if d*d+i+j+2*k+ell-1 == 0 and d*d-a*a+4*k == 0:
                    solutions.append((i,j,k,ell))
check('only identity resonance and no-insertion solution in scan', solutions==[(0,0,0,0),(1,0,0,0)])
# Completeness is analytic: the constant coefficient forces n0+n2+2m+level<=1.

def Asl(a,zz,mu):
    aa=mp.sqrt(1+zz*zz)
    eta=-aa
    x=a-1/(2*zz)
    C=mu/(2*zz*zz)
    return (mp.exp(x/zz*mp.log(C))*mp.gamma(1-2*x*zz)*mp.gamma(1-x/zz)
            *mp.rgamma(mp.mpf('.5')-x*zz+aa*eta)
            *mp.rgamma(mp.mpf('.5')-x*zz-aa*eta))

def RL_source(a,zz,l1):
    # A1 specialization of Fateev2001 eq(25): alpha=sqrt(2), rho=1/sqrt(2).
    bb=-mp.sqrt(2)*zz
    q=(bb+1/bb)/mp.sqrt(2)
    base=l1*mp.gamma(bb*bb)/mp.gamma(1-bb*bb)
    def amp(x):
        u=x-q
        return (mp.exp(u/(mp.sqrt(2)*bb)*mp.log(base))
                *mp.gamma(1-u*mp.sqrt(2)/bb)*mp.gamma(1-u*mp.sqrt(2)*bb))
    return amp(2*q-a)/amp(a)

for tt in [-mp.sqrt(2)/3,-mp.sqrt(3)/4]:
    zz=mp.sqrt(-tt/2)
    aa=zz+1/zz
    mu=mp.mpf('0.63')
    l1=mp.mpf('0.81')
    C=mu/(-tt)
    slope=zz*C**(2/tt-2)*mp.gamma(2-1/tt)/((tt-1)*mp.gamma(1/tt))
    GU=mu/(1-tt)*C**(2/tt-2)*mp.gamma(2-1/tt)/mp.gamma(1/tt)
    GA=(l1*mp.gamma(-tt)/mp.gamma(1+tt))**(1/tt-1)*mp.gamma(2-tt)*mp.gamma(2-1/tt)/(mp.gamma(tt)*mp.gamma(1/tt))
    check('FA Liouville coefficient against rank-one source '+mp.nstr(tt,6), abs(RL_source(-2*zz-1/zz,zz,l1)/GA-1)<mp.mpf('1e-75'))
    for power in [4,6,8]:
        ep=mp.mpf(10)**(-power)
        refl=Asl(1/zz-aa-ep,zz,mu)/Asl(aa+ep,zz,mu)
        err=abs(refl/(ep*slope)-1)
        check('SL reflection has simple zero epsilon=1e-'+str(power)+' t='+mp.nstr(tt,6), err<500*ep)
        check('zero times identity pole gives GU epsilon=1e-'+str(power)+' t='+mp.nstr(tt,6), abs(refl*(-mu/(zz*ep))/GU-1)<500*ep)
    check('candidate GU finite nonzero generic '+mp.nstr(tt,6), mp.isfinite(GU) and abs(GU)>0)
    check('candidate GA finite nonzero generic '+mp.nstr(tt,6), mp.isfinite(GA) and abs(GA)>0)
    s=mp.exp(-2j*mp.pi/tt)
    check('FA coupling-branch monodromy '+mp.nstr(tt,6), abs(mp.exp(2j*mp.pi*(1/tt-1))-1/s)<mp.mpf('1e-75'))
    check('U coupling-branch monodromy '+mp.nstr(tt,6), abs(mp.exp(-2j*mp.pi/tt)-s)<mp.mpf('1e-75'))
    print('DATA t='+mp.nstr(tt,15)+' R_SL/epsilon='+mp.nstr(slope,18)+' GA='+mp.nstr(GA,18)+' GU='+mp.nstr(GU,18))

zero('candidate FA dimensions', (2-2*t)*(1/t-1)-2*(1-t)**2/t)
zero('candidate U dimensions', 1+(2/t-2)-(2/t-1))
zero('candidate AA endpoint dimension', (2-2*t)+2*(1-t)**2/t-(2/t-2))
zero('candidate BB endpoint dimension', 1+(2/t-1)-2/t)
print(f'Completed {count} independent checks.')
