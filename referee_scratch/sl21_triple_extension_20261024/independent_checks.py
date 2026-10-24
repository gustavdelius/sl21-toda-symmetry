"""Fresh local checks for the restricted three-contour and vacuum bridge report.

No earlier check script is imported. These checks supplement the written
distributional argument; they do not construct a quantum operator domain.
"""
import sympy as S
import mpmath as mp

mp.mp.dps = 55
count = 0


def check(name, ok):
    global count
    count += 1
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok:
        raise AssertionError(name)


def eq(name, left, right):
    check(name, S.simplify(S.expand(left - right)) == 0)


def close(name, left, right, tol=mp.mpf("1e-40")):
    check(name, abs(left-right) < tol * max(1, abs(left), abs(right)))


t = S.symbols("t", nonzero=True, real=True)
b = S.sqrt(t)
d0 = S.Matrix([S.I*b/S.sqrt(2), b/S.sqrt(2), S.I])
d2 = S.Matrix([S.I*b/S.sqrt(2), -b/S.sqrt(2), -S.I])
T = d0+d2
d = -T/t
k = d+T
dot = lambda a, b: (a.T*b)[0]
rho = 2/t
beta = t-1
for name, value, target in [
    ("a0 norm", dot(d0,d0), -1),
    ("a2 norm", dot(d2,d2), -1),
    ("mixed screening product", dot(d0,d2), 1-t),
    ("d norm", dot(d,d), -rho),
    ("a0 d", dot(d0,d), 1),
    ("a2 d", dot(d2,d), 1),
    ("A Y2 zero", -dot(-d,d2), 1),
    ("B Y2 exponential pole", -dot(d,d2), -1),
    ("AB exponential power", -dot(-d,d), -rho),
    ("AU endpoint power", -dot(-d,d+d2), 1-rho),
]:
    eq(name, value, target)
for i in range(3):
    eq("Ward derivative vector component %d" % i,
       (d0+beta*d)[i], -(d+d2)[i])

# Rebuild the scalar holomorphic Wick kernel from its two derivative
# contractions, then extract its source poles without using a previous residue.
p, v, r, aa, ww = S.symbols("p v r aa ww", nonzero=True, positive=True)
rr = S.symbols("rho", real=True)
bb = S.symbols("beta", real=True)
hol = (-v)*(v+p)**(-rr)/p*(1/(v+p)+bb/p)
series = S.series(p**2*hol, p, 0, 2).removeO().expand()
V = -v**(1-rr)
eq("source double pole", series.subs(p,0), bb*V)
simple = S.diff(series,p).subs(p,0)
eq("source simple pole is spectator derivative", 
   (simple+S.diff(V,v)).subs({rr:rho,bb:beta}), 0)
# Include the first normal-ordered exponential jet and uncontracted C.
C, D = S.symbols("C D")
full_simple = simple + bb*V*D + C*V
eq("background first jet", full_simple, simple+V*(C+bb*D))
scalar_rp = (p-r)*r**(-rr)/p*(1/r+bb/p)
eq("unexpanded zero-jet decomposition", scalar_rp,
   r**(-rr)*(1/r+(bb-1)/p-bb*r/p**2))

# Integrated superficial degrees derived from the number of independent
# real coordinates and homogeneous factors, after proper Cauchy faces.
degrees = {
    "AB": 1 + (-rr-1),
    "ABY2": 3 + (-rr-1),
    "ABbarB_source_far": 2 + (-rr-1),
    "full_interior": 4 + (-rr-3),
    "full_endpoint": 5 + (-rr-3),
    "tail_cluster": 2 + (-rr-1),
    "AU": 1 + (1-rr),
}
expected = [-rr,2-rr,1-rr,1-rr,2-rr,1-rr,2-rr]
for (name,degree),want in zip(degrees.items(),expected):
    eq("collision degree "+name,degree,want)
    check("degree positive at negative irrational rho "+name,
          bool(degree.subs(rr,-S.sqrt(2))>0))

# The nested first-jet remainder has an exact elementary positive majorant.
# Test the majorant by independent quadrature, including a near-boundary rho.
for rv in [-mp.sqrt(2),-1/mp.sqrt(3)]:
    A, R0 = mp.mpf(5),mp.mpf("0.4")
    for eps in [mp.mpf("1e-3"),mp.mpf("1e-5")]:
        integral = eps*mp.quad(lambda x:x**(-1-rv),[A*eps,R0])
        exact = eps*(R0**(-rv)-(A*eps)**(-rv))/(-rv)
        close("nested remainder quadrature rho=%s eps=%s" % (mp.nstr(rv,6),eps),integral,exact)
        check("nested remainder bounded by constant times epsilon",integral <= eps*R0**(-rv)/(-rv))

# Cochain phases independently generated from the additive labels.
def coc(n,m,tv):
    nL,nR,n0,n2=n
    mL,mR,m0,m2=m
    return (-1)**((nL+nR)*(m0+m2))*mp.exp(-mp.pi*1j*(nL*mR-nR*mL)/tv)

labA=(-1,0,0,0)
labU=(1,1,0,1)
for tv in [-mp.sqrt(2)/3,-mp.sqrt(3)/4]:
    sv=mp.exp(-2j*mp.pi/tv)
    phase=coc(labA,labU,tv)/coc(labU,labA,tv)
    close("AU cocycle ratio",phase,-1/sv)
    close("A right of U exchange",phase*mp.exp(-1j*mp.pi*(1-2/tv)),sv**-2)
    close("A left of U exchange",phase*mp.exp(1j*mp.pi*(1-2/tv)),1)

# These are the coefficients after the separated-face calculation, not a
# proof that arbitrary renormalized products are associative.
s, cc = S.symbols("s cc", nonzero=True)
eq("right endpoint cancels in four-word projection",s**-1*(cc-cc),0)
eq("left endpoint coefficient",s**-1*(-s**2*cc+cc),(s**-1-s)*cc)
eq("raw radial coefficient",t**2-t*(t-1),t)
eq("canonical matched coefficient",t+t*(t-1),t**2)
eq("finite contact propagated in projection",
   (s**-1-s)*(t**2-t),(s**-1-s)*t*(t-1))

# General finite Wick conversion: derive all factors from momenta, with
# the physical two-boson covariance changed and the fermion OPE kept fixed.
c = S.symbols("c", real=True)
H = S.diag(c,c,0)
pU=d+d2
pA=-k
p1=-T
ell = lambda p: S.simplify(-(p.T*H*p)[0]/2) # log Z_p
for name,momentum,want in [
    ("Z0",d0,0),("Z2",d2,0),("Z1",p1,c*t),
    ("ZA",pA,c*(1-t)**2/t),("ZU",pU,c*(1/t-1)),
    ("Zd",d,c/t),
]:
    eq("physical Gaussian "+name,ell(momentum),want)
eq("AA transported endpoint factor",ell(p1)+ell(d)-ell(pA),2*c)
eq("BB transported endpoint factor",ell(d2)+ell(d)-ell(pU),c)
eq("UY0 fusion coefficient",(pU.T*H*d0)[0],c*(1-t))
eq("Y0Y2 fusion coefficient",(d0.T*H*d2)[0],-c*t)
eq("GA covariance at fixed new couplings",ell(pA)+(1/t-1)*ell(p1),c*(1/t-1))
eq("GU covariance at fixed new couplings",ell(pU)+ell(d0)/t+(1/t-1)*ell(d2),c*(1/t-1))
eq("AA fixed-action amplitude scale",2*c-ell(p1)+ell(pA),ell(d))
eq("BB fixed-action amplitude scale",c-ell(d2)+ell(pU),ell(d))

# Direct primary reflection zero, local identity-pole factor, and endpoint
# ratio. Irrational couplings avoid degenerate roots of unity and gamma poles.
def gammaf(x):
    return mp.gamma(x)/mp.gamma(1-x)

for tv in [-mp.sqrt(2)/3,-mp.sqrt(3)/4]:
    zeta=mp.sqrt(-tv/2)
    alpha=mp.sqrt(1+zeta*zeta)
    Q=1/(2*zeta)
    astar=zeta+1/zeta
    l0,l2,l1=mp.mpf("0.73"),mp.mpf("0.91"),mp.mpf("1.17")
    Cstar=mp.sqrt(l0*l2)/(-tv)
    def amp(av):
        v=av-Q
        return Cstar**(v/zeta)*mp.gamma(1-2*v*zeta)*mp.gamma(1-v/zeta)/(mp.gamma(mp.mpf("0.5")-v*zeta-alpha**2)*mp.gamma(mp.mpf("0.5")-v*zeta+alpha**2))
    slope=zeta*Cstar**(2/tv-2)*mp.gamma(2-1/tv)/((tv-1)*mp.gamma(1/tv))
    delta=mp.mpf("1e-15")
    R=amp(2*Q-astar-delta)/amp(astar+delta)
    close("reflection simple-zero slope",R/delta,slope,mp.mpf("1e-12"))
    GU=l0/(1-tv)*Cstar**(2/tv-2)*mp.gamma(2-1/tv)/mp.gamma(1/tv)
    close("zero times identity pole",slope*(-l0/zeta),GU)
    GA=(l1*gammaf(-tv))**(1/tv-1)*mp.gamma(2-tv)*mp.gamma(2-1/tv)/(mp.gamma(tv)*mp.gamma(1/tv))
    AA=(tv/(1-tv))**2*l1*GA
    BB=tv**2*l2*GU
    close("gamma product identity",gammaf(tv)*gammaf(-tv),-1/tv**2)
    close("conditional endpoint ratio",BB/AA,-(-gammaf(tv)*l0*l2/l1)**(1/tv))
    Delta_phi,Delta_theta=1j*mp.pi/zeta,mp.pi/alpha
    directU=mp.exp((zeta+1/zeta)*Delta_phi-1j*alpha*Delta_theta)
    directA=mp.exp(-(2*zeta+1/zeta)*Delta_phi)
    couplingU=mp.exp((1/tv)*(-2j*mp.pi)+(1/tv-1)*0)
    couplingA=mp.exp((1/tv-1)*(2j*mp.pi))
    close("compensating path GU character",couplingU,directU)
    close("compensating path GA character",couplingA,directA)

# Exact translation/coupling covariance, checked by convergent finite
# integrals. This is not a construction of a Toda continuum cycle.
shift=mp.mpf("0.23")*1j
pv=mp.mpf("0.37")
lplus,lminus=mp.mpf("1.2"),mp.mpf("0.8")
def finite_vev(cycle,lp,lm):
    def density(x):
        z=x+cycle
        return mp.exp(-lp*mp.exp(z)-lm*mp.exp(-z))
    den=mp.quad(density,[-12,0,12])
    num=mp.quad(lambda x:mp.exp(pv*(x+cycle))*density(x),[-12,0,12])
    return num/den
direct=finite_vev(shift,lplus,lminus)
pulled=mp.exp(pv*shift)*finite_vev(0,lplus*mp.exp(shift),lminus*mp.exp(-shift))
close("translated cycle uses positive coupling exponents",direct,pulled)
compensated=finite_vev(shift,lplus*mp.exp(-shift),lminus*mp.exp(shift))
close("compensating coupling path has negative exponents",compensated,mp.exp(pv*shift)*finite_vev(0,lplus,lminus))
print("TOTAL %d PASS, 0 FAIL" % count)
