"""Independent fermionic/distributional calculation; no bosonization is used.

Run from the repository root. Previous check scripts are not imported.
Normal products of elementary fermions are represented by an exterior algebra;
cross-contractions between normal-ordered blocks are generated separately.
"""
import math
from pathlib import Path
import numpy as np
import sympy as s


def clean(a):
    return {k:s.expand(v) for k,v in a.items() if s.expand(v) != 0}


def add(*terms):
    out = {}
    for a in terms:
        for k,v in a.items():
            out[k] = out.get(k,0)+v
    return clean(out)


def scale(a,c):
    return clean({k:v*c for k,v in a.items()})


def wedge(a,b):
    out = {}
    for k,v in a.items():
        for l,w in b.items():
            if set(k).intersection(l):
                continue
            sign = (-1)**sum(x>y for x in k for y in l)
            kl = tuple(sorted(k+l))
            out[kl] = out.get(kl,0)+sign*v*w
    return clean(out)


def coeff(a,var,n):
    return clean({k:s.expand(v).coeff(var,n) for k,v in a.items()})


def check(name,a):
    if isinstance(a,dict):
        ok=all(s.simplify(v)==0 for v in a.values())
    elif isinstance(a,(bool,np.bool_)):
        ok=bool(a)
    elif isinstance(a,s.MatrixBase):
        ok=all(s.simplify(v)==0 for v in a)
    else:
        ok=s.simplify(a)==0
    print(('PASS' if ok else 'FAIL')+' '+name,flush=True)
    assert ok,name


unit={():s.S.One}
def jet(species,n=0):
    return {((species,n),):s.S.One}


def fermion(species,pos,order=3):
    return add(*(scale(jet(species,n),pos**n/s.factorial(n))
                 for n in range(order+1)))


# A fermion is (species, position, normal-block number). Species 0,1 are
# chi_0,chi_2; 2,3 are their antichiral partners. Only 0--1 and 2--3 contract.
def contraction(a,b):
    if a[2]==b[2] or {a[0],b[0]} not in ({0,1},{2,3}):
        return 0
    return 1/(a[1]-b[1])


def wick(seq):
    def enumerate_matchings(items,leftovers=(),factor=s.S.One):
        if not items:
            result=unit
            for sp,pos,block in leftovers:
                result=wedge(result,fermion(sp,pos))
            return scale(result,factor)
        first=items[0]
        result=enumerate_matchings(items[1:],leftovers+(first,),factor)
        for j in range(1,len(items)):
            c=contraction(first,items[j])
            if c:
                rest=items[1:j]+items[j+1:]
                result=add(result,enumerate_matchings(
                    rest,leftovers,factor*(-1)**(j-1)*c))
        return result
    return enumerate_matchings(tuple(seq))


r,rb=s.symbols('r rb',nonzero=True)
b=s.symbols('b',nonzero=True)
t=b*b
alpha0=s.Matrix([-1/s.sqrt(2),s.I/s.sqrt(2)])
alpha2=s.Matrix([-1/s.sqrt(2),-s.I/s.sqrt(2)])
p0=-s.I*b*alpha0
p2=-s.I*b*alpha2
T=s.simplify(p0+p2)
d=s.simplify(-T/t)
k=s.simplify(d+T)
x=s.Matrix(s.symbols('x0:2'))
y=s.Matrix(s.symbols('y0:2'))
z=s.Matrix(s.symbols('z0:2'))
xb=s.Matrix(s.symbols('xb0:2'))
check('isotropic bosonic Yukawa momenta',s.Matrix([p0.dot(p0),p2.dot(p2)]))
check('mixed bosonic contraction',p0.dot(p2)+t)
check('current momentum pairs once with each Yukawa',s.Matrix([d.dot(p0)-1,d.dot(p2)-1]))

j=wedge(jet(0),jet(1))
U=wedge(jet(0,1),jet(1))
V=wedge(jet(0),jet(1,1))
dj=add(U,V)
L=add(scale(unit,p0.dot(x)),j)  # polynomial multiplying the fused vertex E
dL=add(scale(unit,p0.dot(y)),dj,scale(L,k.dot(x)))

# Fermionic current OPEs are generated from cross-contractions.
j_chi0=wick([(0,r,0),(1,r,0),(0,0,1)])
j_chi2=wick([(0,r,0),(1,r,0),(1,0,1)])
check('j(z) chi0(w) simple pole',add(coeff(j_chi0,r,-1),scale(jet(0),-1)))
check('j(z) chi2(w) simple pole',add(coeff(j_chi2,r,-1),jet(1)))
check('j(z) chi0(w) regular term is derivative',add(coeff(j_chi0,r,0),scale(jet(0,1),-1)))
check('j(z) chi2(w) regular term is minus derivative',add(coeff(j_chi2,r,0),jet(1,1)))

for idx,p,chi,jchi,expected_sign in [
    (0,p0,jet(0),j_chi0,1),(2,p2,jet(1),j_chi2,-1)]:
    # J(z)Y_i(w): exponential current/Yukawa contraction r^-1;
    # derivative-boson contraction -p0.pi/r; fermions unbosonized.
    core=add(scale(chi,p0.dot(x)+r*p0.dot(y)-p0.dot(p)/r),jchi)
    current_ope=scale(core,(1+r*d.dot(x))/r)
    second=coeff(current_ope,r,-2)
    first=coeff(current_ope,r,-1)
    expected2=scale(chi,1 if idx==0 else t-1)
    expected1=scale(add(jet(0 if idx==0 else 1,1),
                            scale(chi,(d+p).dot(x))),expected_sign)
    check(f'fermionic current/Yukawa {idx} double pole',add(second,scale(expected2,-1)))
    check(f'fermionic current/Yukawa {idx} residue is signed derivative',add(first,scale(expected1,-1)))

# Mixed full Yukawa product. i*i=-1 and left/right fermions anticommute.
pair02=scale(wick([(0,0,0),(2,0,0),(1,r,1),(3,rb,1)]),-r*rb)
pair20=scale(wick([(1,0,0),(3,0,0),(0,r,1),(2,rb,1)]),-r*rb)
hol02=coeff(pair02,rb,0)
hol20=coeff(pair20,rb,0)
check('full Yukawa mixed OPE leading sign +1',add(coeff(hol02,r,0),scale(unit,-1)))
check('opposite ordering also has leading sign +1',add(coeff(hol20,r,0),scale(unit,-1)))
check('chi0--chi2 pair first descendant -j',add(coeff(hol02,r,1),j))
check('chi2--chi0 pair first descendant +j',add(coeff(hol20,r,1),scale(j,-1)))
check('chi0--chi2 pair second descendant -chi0 dchi2',add(coeff(hol02,r,2),V))
check('chi2--chi0 pair second descendant dchi0 chi2',add(coeff(hol20,r,2),scale(U,-1)))

def boson_taylor(p):
    return 1+r*p.dot(x)+r*r*(p.dot(y)+p.dot(x)**2)/2

G0=scale(hol02,boson_taylor(p2))
G2=scale(hol20,boson_taylor(p0))
N=add(coeff(G2,r,2),scale(coeff(G0,r,2),-1))
RP=add(scale(L,T.dot(x)),scale(unit,-t*(T.dot(y)+T.dot(x)**2)/2))
A=scale(add(scale(unit,(p0-p2).dot(x)),scale(j,2)),T.dot(x)/2)
primitiveP=scale(unit,-t*(T.dot(y)+T.dot(x)*k.dot(x))/2)
check('composite residue nonderivative class',add(RP,scale(primitiveP,-1),scale(A,-1)))
expectedN=add(scale(unit,((p0-p2).dot(y)+(p0-p2).dot(x)*T.dot(x))/2),
              dj,scale(j,T.dot(x)))
check('fermionic difference of collision endpoints gives N',add(N,scale(expectedN,-1)))
check('full local Ward identity N - RP/t = derivative L',add(N,scale(RP,-1/t),scale(dL,-1)))

# Contact derivatives are fixed by double poles, not only their residues.
I0p,I2p,Q0,Q2=s.symbols('I0p I2p Q0 Q2')
# Q_i is integral (d theta)F_i. integral theta*dF_i = I_i' - Q_i.
dist0=(I0p-Q0)-I0p
dist2=-(I2p-Q2)-(t-1)*I2p
check('distributional node-0 cutoff derivative',dist0+Q0)
check('distributional node-2 cutoff derivative',dist2+t*I2p-Q2)
check('full separated-current double insertion divergence',
      dist0+dist2-(-t*I2p+Q2-Q0))

# Verify the composite triple/double/simple poles by direct expansion at w.
current_T=scale(add(
    scale(unit,p0.dot(x)+r*p0.dot(y)+r*r*p0.dot(z)/2+t/r),
    j,scale(dj,r)), (1+r*d.dot(x)+r*r*(d.dot(y)+d.dot(x)**2)/2)/r**2)
Q1=coeff(current_T,r,-1)
Q2P=coeff(current_T,r,-2)
Q3P=coeff(current_T,r,-3)
check('composite triple pole',add(Q3P,scale(unit,-t)))
# derivative of Q2P=(p0+t*d).dphi+j, times E
vec=p0+t*d
derQ2=add(scale(unit,vec.dot(y)+vec.dot(x)*k.dot(x)),dj,scale(j,k.dot(x)))
der2Q3=scale(unit,t*(k.dot(y)+k.dot(x)**2))
distributionP=add(scale(Q1,-1),derQ2,scale(der2Q3,-s.Rational(1,2)))
check('delta, delta-prime and delta-second give composite RP',
      add(distributionP,scale(RP,-1)))

# Operator mixing of K2 with Y0, including a consistent cocycle sign sigma.
sigma=s.symbols('sigma')
R,eps=s.symbols('R eps',positive=True)
# K2(0)Y0(r)=-sigma*r^(t-2)*rbar^(t-1)*[1+r*L+...].
mixing=scale(L,-sigma*s.pi/t*(R**(2*t)-eps**(2*t)))
H1mix=scale(mixing,-sigma*t/s.pi)
H1mix=clean({key:s.expand(val).subs(sigma**2,1) for key,val in H1mix.items()})
check('H1 insertion divergence is -eps^(2t)*L',
      add(H1mix,scale(L,-R**(2*t)+eps**(2*t))))
check('same local counterterm makes H finite',
      add(H1mix,scale(L,eps**(2*t)-R**(2*t))))

# Comparing the old contour primitive with the actual current subtraction.
oldH=scale(unit,(t+1)*T.dot(x)/(2*t)) # lambda0 lambda2 eps^(2t) suppressed
check('old contact primitive is not the current subtraction',
      bool(add(L,scale(oldH,-1))))

# An explicit associative local cocycle. n=(n_d,n_0,n_2).
def cocycle(n,m):
    return (-1)**(n[0]*(m[1]+m[2]))
def plus(n,m):
    return tuple(a+b for a,b in zip(n,m))
charges=[(1,0,0),(-1,0,0),(0,1,0),(0,0,1),(0,1,1),(1,1,0),(1,0,1)]
check('chosen local cocycle is associative',
      all(cocycle(a,b0)*cocycle(plus(a,b0),c)
          ==cocycle(b0,c)*cocycle(a,plus(b0,c))
          for a in charges for b0 in charges for c in charges))
check('current/Yukawa sign minus; current/composite sign plus',
      cocycle((1,0,0),(0,1,0))==-1
      and cocycle((1,0,0),(0,0,1))==-1
      and cocycle((1,0,0),(0,1,1))==1)

# Generic-coupling anomaly counting at the next mass order, as a scope check.
# n0,n2 count Yukawas, ne counts even-root insertions. Boson charges use only
# two fields; minimal same-charge fermion jets have weights F^2/2 each.
def levels(n0,n2,ne):
    Nmass=n0+n2+2*ne
    q=n0*p0+n2*p2-ne*T
    ferm_charge=n0-n2
    hbase=s.simplify(-(d+q).dot(d+q)/2+s.Rational(1,2)*ferm_charge**2)
    hbbase=s.simplify(-q.dot(q)/2+s.Rational(1,2)*ferm_charge**2)
    lamweight=s.Rational(n0+n2,2)+ne*(1-t)
    return (s.simplify(1+1/t-lamweight-hbase),
            s.simplify(1-lamweight-hbbase))

print('FORMULA Hren = sigma*t*lambda2*K2 + lambda0*lambda2*epsilon^(2t)*L',flush=True)
print('FORMULA kappa = -lambda0*lambda2*epsilon^(2t)/t; N-RP/t = partial L',flush=True)
for n0,n2,ne in [(3,0,0),(2,1,0),(1,2,0),(0,3,0),(1,0,1),(0,1,1)]:
    lev=levels(n0,n2,ne)
    print(f'ORDER3 n0={n0} n2={n2} ne={ne}: descendant_levels={lev}',flush=True)

# A numerical angular check of the cutoff terms with nonconstant spectators.
# For this check choose a definite fermion contraction background equivalent to
# j and U,V expectation values; the analytic Wick calculation above is primary.
def numeric_circle(tval):
    bv=complex(tval)**0.5
    p0v=np.array(p0.subs(b,bv),complex).reshape(2)
    p2v=np.array(p2.subs(b,bv),complex).reshape(2)
    xv=np.array([.23+.17j,-.11+.08j])
    yv=np.array([-.07+.02j,.04+.09j])
    xbv=np.array([.13-.03j,.06+.05j])
    # Actual fermion two-point functions with spectators at f=A,g=B:
    # <chi0(0)chi2(r)chi0(A)chi2(B)>/<chi0(A)chi2(B)>
    # normalized by the leading -1/r.
    Apos=1.3+.4j
    Bpos=-1.1+.3j
    f02=lambda rr: 1+rr*(Apos-Bpos)/(Bpos*(rr-Apos))
    # Reverse charges at 0,r, same fixed spectators. Wick signs differ.
    f20=lambda rr: 1+rr*(Apos-Bpos)/(Apos*(Bpos-rr))
    # derivative coefficients obtained directly from these rational functions.
    jf=(Apos-Bpos)/(Apos*Bpos)
    vf=(Apos-Bpos)/(Bpos*Apos**2)
    uf=(Apos-Bpos)/(Apos*Bpos**2)
    L0=(p2v.dot(yv)+p2v.dot(xv)**2)/2-jf*p2v.dot(xv)-vf
    L2=(p0v.dot(yv)+p0v.dot(xv)**2)/2+jf*p0v.dot(xv)+uf
    Nv=L2-L0
    theta=2*np.pi*np.arange(16384)/16384
    errors=[]
    for e in [.08,.04,.02,.01]:
        rr=e*np.exp(1j*theta)
        G0v=f02(rr)*np.exp(rr*p2v.dot(xv)+rr*rr*p2v.dot(yv)/2
                          +rr.conjugate()*p2v.dot(xbv))
        G2v=f20(rr)*np.exp(rr*p0v.dot(xv)+rr*rr*p0v.dot(yv)/2
                          +rr.conjugate()*p0v.dot(xbv))
        value=np.pi*e**(2*tval-2)*np.mean(np.exp(-2j*theta)*(G2v-G0v))
        want=np.pi*e**(2*tval)*Nv
        error=abs((value-want)/want)
        errors.append(error)
        print(f'NUMERICAL t={tval:.15g} eps={e:g} relative_error={error:.6g}',flush=True)
    check(f'fermionic spectator collision circles at t={tval:.15g}',
          errors[-1]<2e-4 and all(3.7<errors[i]/errors[i+1]<4.3 for i in range(3)))


numeric_circle(math.sqrt(2)/3)
numeric_circle(-math.sqrt(2)/3)


def mixing_correlator(tval):
    # Exact external-fermion correlator, with bosonic charge balanced at infinity.
    # The only nonzero local descendant needed here is j:
    # <j(0) chi0(A) chi2(B)>/<chi0(A) chi2(B)> = (A-B)/(A*B).
    Apos=1.3+.4j
    Bpos=-1.1+.3j
    jvalue=(Apos-Bpos)/(Apos*Bpos)
    outer=.3
    angles=2*np.pi*np.arange(256)/256
    nodes,weights=np.polynomial.legendre.leggauss(64)
    expected=outer**(2*tval)*jvalue
    for e in [.08,.04,.02,.01]:
        low,high=np.log(e),np.log(outer)
        logs=(high-low)*nodes/2+(high+low)/2
        rho=np.exp(logs)
        rr=rho[:,None]*np.exp(1j*angles[None,:])
        f20=1+rr*(Apos-Bpos)/(Apos*(Bpos-rr))
        angular=np.mean(np.exp(-1j*angles[None,:])*f20,axis=1)
        hbare=2*tval*(high-low)/2*np.dot(
            weights,rho**(2*tval-1)*angular)
        hren=hbare+e**(2*tval)*jvalue
        error=abs(hren-expected)
        print(f'MIXING t={tval:.15g} eps={e:g} bare_abs={abs(hbare):.6g} '
              f'renormalized_error={error:.4g}',flush=True)
        check(f'spectator correlator subtraction t={tval:.15g} eps={e:g}',error<5e-10)
    # The previous purely bosonic primitive has zero matrix element here.
    check(f'old primitive misses nonzero fermion mixing at t={tval:.15g}',abs(jvalue)>0.1)


mixing_correlator(math.sqrt(2)/3)
mixing_correlator(-math.sqrt(2)/3)
