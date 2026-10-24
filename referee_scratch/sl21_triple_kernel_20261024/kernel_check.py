"""Independent first-Yukawa three-current Wick and contour checks.

No project or earlier referee script is imported. Lambda_2 is set to one.
The exact unexpanded numerical kernels use A(a) B(z) bar B(0), a=3.
The ordering cocycle makes the common A-U times B-barB phase equal to -1.
"""

import cmath
import math
import sympy as S
from scipy.integrate import quad

checks = 0


def check(label, condition, detail=""):
    global checks
    checks += 1
    print(f"{'PASS' if condition else 'FAIL'} {label}" + (f": {detail}" if detail else ""), flush=True)
    if not condition:
        raise AssertionError(label)


def zero(expr):
    return S.simplify(S.factor(S.together(expr))) == 0


t, r, p, a, z, w = S.symbols("t r p a z w", nonzero=True)
rho, beta = 2/t, t-1
scalar = (p-r)/p * (1/r + beta/p)
check("exact zero-jet chiral Wick polynomial", zero(scalar - (1/r+(beta-1)/p-beta*r/p**2)))

# This contraction is a direct consequence of the bilinear root products.
G = (a-z)**(-rho)*(a-w)
double = S.simplify(beta*G.subs(z,w))
simple = S.simplify((G/(z-a)+beta*S.diff(G,z)).subs(z,w))
V = (a-w)**(1-rho)
check("mixed face double pole is beta V", zero(double-beta*V))
check("mixed face simple residue is minus V prime", zero(simple+S.diff(V,w)))

D, C, Aleft = S.symbols("D C Aleft")
# The physical root identity is C+(t-1)D=-(D+Aleft).
Croot = -t*D-Aleft
Deff = D+rho/a
Ceff = Croot-1/a
Aeff = Aleft-1/a
Keff = Deff+Aeff
check("separated A preserves the full first-jet Ward identity", zero(Ceff+beta*Deff+Keff))
check("dressed endpoint derivative", zero(S.diff(V,w).subs(w,0)/V.subs(w,0)-(rho-1)/a))

# Exact coefficients of the two tail singularities, before the conjugate
# Plemelj signs. Their common cocycle is -1 relative to the BBbar formula.
tail1 = -t*a*(a-z)**(-rho)*(1/((z-a)*z)+beta/z**2)
V0 = a**(1-rho)
V1 = (rho-1)*a**(-rho)
tail1_double = S.limit(z**2*tail1,z,0)
tail1_simple = S.diff(S.cancel(z**2*tail1),z).subs(z,0)
check("first tail double coefficient", zero(tail1_double+t*beta*V0))
check("first tail simple coefficient", zero(tail1_simple-t*V1))
check("tail contour coefficient is t squared", zero(t+t*beta-t**2))
check("action plus tails gives radial t", zero(t**2-t*beta-t))
check("finite matching restores t squared", zero(t+t*beta-t**2))

# Genuine irrational coupling, with a common first-order convergence region.
tn = -math.sqrt(2)
rn = 2/tn
bn = tn-1
an = 3.0
R = 1.0
cn = math.sqrt(3)
v = lambda ww: (an-ww)**(1-rn)
vp = lambda ww: -(1-rn)*(an-ww)**(-rn)
v1 = vp(0.0)


def cq(fun, lo, hi, points=None, tol=2e-9):
    kw = dict(epsabs=tol, epsrel=tol, limit=250)
    if points:
        kw['points'] = [q for q in points if lo < q < hi]
    return quad(lambda xx: complex(fun(xx)).real, lo, hi, **kw)[0] + 1j*quad(lambda xx: complex(fun(xx)).imag, lo, hi, **kw)[0]


def action_integrals(eps):
    delta=cn*eps
    # Integrate the unexpanded V'(w), with an actual circular inner hole.
    # Splitting the y integral at +-eps retains the moving hole boundary.
    def horizontal(y):
        outer=math.sqrt(R*R-y*y)
        inner=math.sqrt(max(0.0,eps*eps-y*y))
        fn=lambda x: vp(complex(x,y))/complex(x,-y)**2
        if inner:
            return cq(fn,-outer,-inner)+cq(fn,inner,outer)
        return cq(fn,-outer,outer,points=[0.0])
    I=cq(horizontal,-delta,delta,points=[-eps,0.0,eps],tol=2e-8)
    edge=math.sqrt(R*R-delta*delta)
    J=(-0.5j)*cq(lambda x: v(complex(x,-delta))/complex(x,delta)**2-v(complex(x,delta))/complex(x,-delta)**2,-edge,edge,points=[0.0])
    # All insertion-order cocycles and the action coefficient are included.
    return 2j*bn*(I-bn*J)


def tail_integral(eps):
    delta=cn*eps
    def upper(x):
        zz=complex(x,delta)
        return -tn*an*(an-zz)**(-rn)*(1/((zz-an)*zz)+bn/zz**2)-tn*bn*v(zz)/zz.conjugate()**2
    return 2j*cq(lambda x: upper(x).imag,-R,R,points=[0.0])


target_action=2j*math.pi*tn*bn*v1
target_tail=-2j*math.pi*tn**2*v1
target_radial=-2j*math.pi*tn*v1
rows=[]
for eps in [0.02,0.01,0.005,0.0025,0.00125,0.000625]:
    action=action_integrals(eps)
    tails=tail_integral(eps)
    radial=action+tails
    rows.append((eps,action,tails,radial))
    print(f"DATA eps={eps:.5f} action={action.imag:.12f} tails={tails.imag:.12f} radial={radial.imag:.12f}",flush=True)
print(f"DATA targets action={target_action.imag:.12f} tails={target_tail.imag:.12f} radial={target_radial.imag:.12f}",flush=True)
for col,label,target in [(1,'action',target_action),(2,'tails',target_tail),(3,'radial sum',target_radial)]:
    errors=[abs(row[col]-target) for row in rows]
    check(f"unexpanded {label} contour converges",all(errors[k+1]<0.55*errors[k] for k in range(len(rows)-1)),f"last error {errors[-1]:.6g}")
    extrap=2*rows[-1][col]-rows[-2][col]
    check(f"unexpanded {label} extrapolation",abs(extrap-target)<2e-4,f"error {abs(extrap-target):.6g}")


# Independent numerical scaling of the full scalar kernel, away from proper
# faces. Multiplication by eps^4 is the relative interior measure.
def raw_kernel(aa,zz,uu,ww):
    return (aa-zz)**(-rn)*(aa-ww)/(zz-ww)*(1/(zz-aa)+bn/(zz-ww))*bn/(uu.conjugate()-ww.conjugate())**2
coords=(complex(-1.7,0.4),complex(0.8,0.9),complex(0.2,-0.4),complex(-0.3,1.8))
base=raw_kernel(*coords)
for scale in [0.5,0.25,0.125]:
    val=raw_kernel(*(scale*q for q in coords))*scale**4/base
    expected=scale**(1-rn)
    check(f"interior simultaneous homogeneity at scale {scale}",abs(val-expected)<1e-11,f"factor={val.real:.12g}, degree={1-rn:.12g}")
check("proper AB face converges",-rn>0)
check("simultaneous interior face converges",1-rn>0)
check("simultaneous endpoint face converges",2-rn>0)
check("A with fused endpoint converges",2-rn>0)

print(f"TOTAL {checks} PASS, 0 FAIL",flush=True)
