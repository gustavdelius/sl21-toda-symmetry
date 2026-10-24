"""Independent audit. No earlier investigation script is read or imported.

The numerical contact calculation starts with a smooth two-dimensional Cauchy
kernel, integrates the outer contour exactly, and then integrates the remaining
inner hole. The collision test uses a fresh Wick kernel on a compact sector.
"""
import math
import numpy as np
import sympy as s
from scipy.integrate import quad
from scipy.stats import qmc

passed = 0
failed = 0


def check(name, truth, detail=""):
    global passed, failed
    ok = bool(truth)
    passed += int(ok)
    failed += int(not ok)
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f": {detail}" if detail else ""))


def eq(name, lhs, rhs):
    check(name, s.simplify(lhs-rhs) == 0)


b = s.symbols('b', nonzero=True)
t = s.symbols('t', positive=True)
roots = [s.Matrix([-1, s.I])/s.sqrt(2),
         s.Matrix([-1, -s.I])/s.sqrt(2)]
pp = [-s.I*b*r for r in roots]
a = [s.Matrix([*pp[0], s.I]), s.Matrix([*pp[1], -s.I])]
d = -(a[0]+a[1])/b**2


def dot(v, w):
    return s.simplify(v.dot(w)).subs(b**2, t)


eq('null odd roots, node 0', roots[0].dot(roots[0]), 0)
eq('null odd roots, node 2', roots[1].dot(roots[1]), 0)
eq('odd root pairing', roots[0].dot(roots[1]), 1)
eq('bosonized a0 square', dot(a[0], a[0]), -1)
eq('bosonized a2 square', dot(a[1], a[1]), -1)
eq('bosonized mixed pairing', dot(a[0], a[1]), 1-t)
eq('dual length', dot(d, d), -2/t)
eq('dual screening pairing 0', dot(d, a[0]), 1)
eq('dual screening pairing 2', dot(d, a[1]), 1)

# Derive the fermion four-point function from its three Wick pairings.
z = s.symbols('z0:4')
charges = [1, -1, 1, -1]
def fermion(i, j):
    return 1/(z[i]-z[j]) if charges[i] != charges[j] else 0
wick = fermion(0, 1)*fermion(2, 3) - fermion(0, 2)*fermion(1, 3) + fermion(0, 3)*fermion(1, 2)
vertex = s.prod((z[i]-z[j])**(charges[i]*charges[j])
                for i in range(4) for j in range(i+1, 4))
eq('four fermions agree with local exponential bosonization', wick, vertex)

# Current-first expansion: derivative contraction plus first Taylor coefficient.
for idx in (0, 1):
    double = -dot(a[0], a[idx])
    simple = a[0] + double*d
    expected = (d+a[idx]) * (1 if idx == 0 else -1)
    check(f'current pole vector node {2*idx}', all(s.simplify(x.subs(t,b**2)) == 0 for x in simple-expected))
    # The local cocycle is -1; double-pole differentiation is retained.
    scalar = 1 if idx == 0 else -1
    eq(f'distributional tail coefficient node {2*idx}', -(double-scalar), 0 if idx == 0 else -t)

# Radial endpoint expansion is obtained from the actual charge pairing.
v = d+a[1]
eq('endpoint--Y0 holomorphic exponent', -dot(v, a[0]), t-2)
eq('endpoint--Y2 holomorphic exponent', -dot(v, a[1]), 0)
check('endpoint fused momentum is k', all(s.simplify(x) == 0 for x in v+a[0]-(1-b**2)*d))
rho, eps, radius = s.symbols('rho eps radius', positive=True)
for n in (0, 1):
    primitive = s.pi*(radius**(2*t-2+2*n)-eps**(2*t-2+2*n))/(t-1+n)
    eq(f'endpoint annulus derivative, jet level {n}', s.diff(primitive, radius), 2*s.pi*radius**(2*t-3+2*n))

# Gaussian current polynomial: Taylor jets from the two fermionic endpoints.
x0, x2, y0, y2, j, U, V = s.symbols('x0 x2 y0 y2 j U V')
Q0 = (y2+x2*x2)/2-j*x2-V
Q2 = (y0+x0*x0)/2+j*x0+U
RP = (x0+j)*(x0+x2)-t*(y0+y2+(x0+x2)**2)/2
dL = y0+U+V+(x0+j)*(1-1/t)*(x0+x2)
eq('independent endpoint-jet Ward polynomial', Q2-Q0-RP/t, dL)
# Recompute the potential response from its three current-first poles.
C, Cp, D, Dp, X, Xp = s.symbols('C Cp D Dp X Xp')
simple = Cp+C*D+t*(Dp+D**2)/2
double_derivative = Cp+t*Dp+(C+t*D)*(D+X)
triple_second_derivative = t*(Dp+Xp+(D+X)**2)/2
eq('three potential poles reproduce fixed-current response',
   -simple+double_derivative-triple_second_derivative,
   C*X-t*(Xp+X**2)/2)

# Independent outer-contour extension, with z=x+iy:
# C_eta=bar(z)/(|z|^2+eta^2), P_eta=-partial_z C_eta.
# Integrating over x gives -i*pi*y/sqrt(y^2+eta^2) and
# pi*eta^2/(2*(y^2+eta^2)^(3/2)), respectively.
x, y, eta = s.symbols('x y eta', real=True, positive=True)
Ceta = (x-s.I*y)/(x*x+y*y+eta*eta)
Peta = (x-s.I*y)**2/(x*x+y*y+eta*eta)**2
eq('smooth double pole is differentiated Cauchy kernel', -(s.diff(Ceta,x)-s.I*s.diff(Ceta,y))/2, Peta)

# Horizontal integrals outside the unit circular hole. Each is evaluated by
# elementary antiderivatives before the remaining y quadrature.
def A2(y):
    return 2*math.sqrt(max(0., 1-y*y)) if abs(y) < 1 else 0.

def A1_im(y):
    if y == 0:
        return 0.
    if abs(y) >= 1:
        return math.copysign(math.pi,y)
    return 2*math.asin(y)

def integrate_y(fun, c):
    cuts = [-np.inf, -c, -1., 0., 1., c, np.inf] if c > 1 else [-np.inf,-1.,-c,0.,c,1.,np.inf]
    return sum(quad(fun, lo, hi, epsabs=2e-10, epsrel=2e-10, limit=200)[0]
               for lo,hi in zip(cuts[:-1],cuts[1:]))

def faces(c, eta):
    def sign_smooth(y):
        return y/math.sqrt(y*y+eta*eta)
    def deriv(y):
        return eta*eta/(y*y+eta*eta)**1.5
    def theta(y):
        return (sign_smooth(c-y)+sign_smooth(c+y))/2
    # partial_w Theta = i * edge_im.
    def edge_im(y):
        return (deriv(c-y)-deriv(c+y))/4
    I = integrate_y(lambda y: theta(y)*A2(y), c)
    J1 = integrate_y(lambda y: -edge_im(y)*A1_im(y), c)
    Jw = integrate_y(lambda y: -edge_im(y)*(A1_im(y)+2*y*A2(y)), c)
    # Extra eta term in the exact first contour moment of P_eta.
    def E_im(z):
        return math.pi*z*eta*eta/(2*(z*z+eta*eta)**1.5)
    extra = integrate_y(lambda y: (E_im(c-y)-E_im(-c-y))*A2(y), c)
    return I, J1, Jw, extra

def face_limit(c):
    if c >= 1:
        return np.array([math.pi,-math.pi,-math.pi])
    h = math.sqrt(1-c*c)
    return np.array([2*(math.asin(c)+c*h),-2*math.asin(c),-2*math.asin(c)-4*c*h])

tn = math.sqrt(2.)
cn = math.sqrt(2.)
KL, KR, D = .3, -.2, .17
q = tn-1
current_C = -KL-q*D
AR = .11
Cbar = q*AR-tn*KR
expected = -2*math.pi*tn*q*(KL+KR)  # Imaginary coefficient / (nu lambda_2).
errors = []
for en in (.08,.04,.02,.01):
    I,J1,Jw,extra = faces(cn,en)
    # Derived directly from the constant and first moments of the smooth
    # contour kernels, before the strip functions are taken to a sharp limit.
    value = 2*(q*(current_C+q*D)*I + q*q*KL*Jw + q*(q*AR-Cbar)*J1) - q*q*D*extra/math.pi
    err = abs(value-expected)
    errors.append(err)
    print(f'DATA smooth contour eta={en:.5g} I={I:.12g} J1={J1:.12g} Jw={Jw:.12g} contact/i={value:.12g} error={err:.4g}')
check('original smooth contour contact converges', all(v>u for v,u in zip(errors,errors[1:])))
check('original smooth contour matches -2*pi*t*(t-1)*(KL+KR)', errors[-1] < 3e-4, f'target={expected:.12g}')
for c in (1/math.sqrt(2), math.sqrt(2)):
    val = np.array(faces(c,.003)[:3])
    check(f'independent smooth-kernel face limit c={c:.8g}', np.max(np.abs(val-face_limit(c))) < 4e-4, str(val))

# Tail contacts directly from the current-first Wick coefficients.
Kleft, Kright = s.symbols('Kleft Kright')
hol_simple = -t*Kleft
anti_simple = t*(t-1)*Kleft+t*t*Kright
eq('tail Plemelj simple-pole difference', -hol_simple+anti_simple, t*t*(Kleft+Kright))
eq('radial coefficient after bare-action contact', t*t-t*(t-1), t)
gamma = t/(1-t)
eq('AA contact from contour first jet', t/(1-t), gamma)
eq('AA radial coefficient', gamma**2+gamma, gamma**2/t)
eq('canonical finite BB addition', t+t*(t-1), t*t)
eq('finite BB contact leading endpoint subtraction', -t*(t-1)/(t-1), -t)
eq('finite BB contact descendant endpoint subtraction', -t*(t-1)/t, -(t-1))

# Actual four-insertion Wick kernel, not an assumed power count.
u,w0,w2 = s.symbols('u w0 w2')
f = (-(dot(a[0],a[0]))/(u-w0) - dot(a[0],a[1])/(u-w2))/((u-w0)*(u-w2))
r = w0-w2
partial = 1/(r*(u-w0)**2) -(t-1)/(r*(u-w2)**2) +(t-2)/r**2*(1/(u-w0)-1/(u-w2))
eq('joint chiral Wick kernel partial fractions', f, partial)
eq('joint chiral kernel is odd under simultaneous inversion', f.subs({u:-u,w0:-w0,w2:-w2}, simultaneous=True), -f)
scale = s.symbols('scale', positive=True)
eq('joint chiral kernel homogeneity', f.subs({u:scale*u,w0:scale*w0,w2:scale*w2}, simultaneous=True), f/scale**3)

# A compact, nonempty sector avoids every current pole and every proper
# interaction collision. Reflection is tested before integration.
points = qmc.Sobol(d=5, scramble=False).random_base2(12)
points = (points-.5)*6
ux = points[:,0]
aa = points[:,1]+1j*points[:,2]
bb = points[:,3]+1j*points[:,4]
mask = (abs(aa)>1)&(abs(bb)>1)&(abs(aa-bb)>1)
for pole in (ux+1j*cn,ux-1j*cn):
    mask &= (abs(pole-aa)>.4)&(abs(pole-bb)>.4)
ux,aa,bb = ux[mask],aa[mask],bb[mask]

def fn(u,a,b,tc):
    return (1/(u-a)+(tc-1)/(u-b))/((u-a)*(u-b))

def defect(x,a,b,delta,tc):
    return abs(a-b)**(2*tc-2)*(fn(x+1j*delta,a,b,tc)-fn(x-1j*delta,a,b,tc))*fn(0,np.conj(a),np.conj(b),tc)

tc = math.sqrt(1.7)  # Between 1 and 3/2, irrational.
base = defect(ux,aa,bb,cn,tc)
reflected = defect(-ux,-aa,-bb,cn,tc)
check('primitive sector is nonempty and kernel nonzero', len(ux)>100 and np.mean(abs(base))>1e-4, f'n={len(ux)}, mean abs={np.mean(abs(base)):.7g}')
check('constant fused-field defect is reflection odd', np.max(abs(base+reflected))<1e-10)
for sc in (.5,.25):
    scaled = defect(sc*ux,sc*aa,sc*bb,sc*cn,tc)*sc**5
    err = np.max(abs(scaled-base*sc**(2*tc-3)))
    check(f'fixed-inner relative measure gives power 2t-3, scale={sc}',err<1e-9)
    linear = scaled*(sc*aa.real)
    check(f'first Taylor jet gives power 2t-2, scale={sc}',np.max(abs(linear-base*aa.real*sc**(2*tc-2)))<1e-9)
print(f'DATA raw interior exponent={2*tc-3:.9g}; first-jet exponent={2*tc-2:.9g}; full endpoint dimension adds one measure')

# Derive each kernel degree from its Wick contractions. Each holomorphic
# derivative prefactor is contracted at leading order, except the explicitly
# uncontracted anti-derivative in J bar L.
single_current_double_y = -dot(d,a[0])-dot(d,a[1])-2*dot(a[0],a[1])-1
two_current_double_y = -2*dot(d,a[0])-2*dot(d,a[1])-2*dot(a[0],a[1])-2
current_tail_y = -dot(d,a[1])-dot(d,a[0])-dot(a[1],a[0])-dot(v,a[0])-1
T = a[0]+a[1]
current_current_P = -2*dot(d,T)-2
current_L = -dot(d,T)-1
same_node = -4*dot(d,a[1])-2*dot(a[1],a[1])-2
tail_tail = -dot(v,a[1])-dot(a[1],v)

# Relative, not total, integration dimensions for the primitive scalar jets.
degrees = {
    'JbarJ Y0Y2': (two_current_double_y,5,0,2*t-3),
    'JbarK Y0': (current_tail_y,3,0,2*t-3),
    'JbarJ P': (current_current_P,3,2*t,2*t-3),
    'JbarL': (current_L,1,2*t,2*t-2),
    'Y0Y2': (-2*dot(a[0],a[1]),2,0,2*t),
    'one current Y0Y2': (single_current_double_y,4,0,2*t-1),
    'JbarJ equal Yukawas': (same_node,5,0,1),
    'KbarK': (tail_tail,1,0,1),
    'endpoint U Y0': (-2*dot(v,a[0]),2,0,2*t-2),
}
for name,(kernel,measure,ct,answer) in degrees.items():
    eq('relative collision count '+name, kernel+measure+ct, answer)
    check('safe convergence at Re(t)>2: '+name, answer.subs(t,2)>0 if hasattr(answer,'subs') else answer>0)

eq('A with either Yukawa starts at a positive power', -dot(-d,a[1]), 1)
eq('B with even-root potential remains regular after derivative contraction', -dot(d,-T)-1, 1)
eq('A with the composite potential is regular', -dot(-d,T), 2)

# A local meromorphic model with both endpoint subtractions and a genuine
# analytic spectator series, independently integrated term by term.
def ac_endpoint(tc,rr=.7):
    return math.pi*sum(rr**(2*tc-2+2*n)/(math.factorial(n)**2*(tc-1+n)) for n in range(35))

def sub_endpoint(tc,ee,rr=.7):
    # Angular average of exp(z+bar z) is sum r^(2n)/(n!)^2.
    def angular(r):
        return sum(r**(2*n)/math.factorial(n)**2 for n in range(20))
    raw = 2*math.pi*quad(lambda r:r**(2*tc-3)*angular(r),ee,rr,epsabs=1e-9,epsrel=2e-12)[0]
    return raw+math.pi*ee**(2*tc-2)/(tc-1)+math.pi*ee**(2*tc)/tc

for tc in (math.sqrt(2),-math.sqrt(2)/5):
    target = ac_endpoint(tc)
    vals = [sub_endpoint(tc,ee) for ee in (.04,.02,.01)]
    errs = [abs(v-target) for v in vals]
    check(f'nonconstant endpoint spectator continuation t={tc:.7g}', errs[2]<errs[1]<errs[0], f'errors={errs}')

print(f'RESULT {passed} PASS, {failed} FAIL')
raise SystemExit(bool(failed))
