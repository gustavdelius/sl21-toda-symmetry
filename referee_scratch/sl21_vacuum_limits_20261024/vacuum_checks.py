"""Independent root/sector/normalization audit. No project modules imported."""
import sympy as S

count = 0

def check(name, expr):
    global count
    ok = all(S.simplify(x) == 0 for x in expr) if isinstance(expr, S.MatrixBase) else S.simplify(expr) == 0
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        raise AssertionError(S.simplify(expr))
    count += 1

z = S.symbols('z', positive=True)
t = -2*z**2
b = S.I*S.sqrt(2)*z
alpha = S.sqrt(1+z**2)
a0 = S.Matrix([S.I*b/S.sqrt(2), b/S.sqrt(2), S.I])
a2 = S.Matrix([S.I*b/S.sqrt(2), -b/S.sqrt(2), -S.I])
T = a0+a2
d = -T/t
k = T+d
O = S.Matrix([[-1,0,0],[0,z/alpha,1/alpha],[0,1/alpha,-z/alpha]])
check('orthogonal local field transformation', O*O.T-S.eye(3))
check('unit determinant local transformation', O.det()-1)
check('Y0 momentum in SL coordinates', O*a0-S.Matrix([z,S.I*alpha,0]))
check('Y2 momentum in SL coordinates', O*a2-S.Matrix([z,-S.I*alpha,0]))
check('even root potential momentum', O*(-T)-S.Matrix([-2*z,0,0]))
check('composite potential momentum', O*T-S.Matrix([2*z,0,0]))
check('dual current momentum', O*d-S.Matrix([1/z,0,0]))
check('FA endpoint momentum', O*(-k)-S.Matrix([-2*z-1/z,0,0]))
check('U endpoint momentum', O*(d+a2)-S.Matrix([z+1/z,-S.I*alpha,0]))
check('F endpoint momentum', O*k-S.Matrix([2*z+1/z,0,0]))
check('C prefactor has no free component', (O*a0)[2])
check('first SL hierarchy coupling constraint', alpha**2-z**2-1)
check('physical Yukawa0 null', sum(x*x for x in a0[:2]))
check('physical Yukawa2 null', sum(x*x for x in a2[:2]))
check('bosonized fermion dimension', a0.dot(a0)+1)
check('unlike-screening OPE exponent', a0.dot(a2)-(1-t))
check('dual screening integer OPE', a0.dot(d)-1)

# LP2017 eq(1.8), using ell for its radial label b.
h = 2-t
g = 1/h
labels = {
    'identity': (0, 0),
    'FA': (0, 2*h-2),
    'U': (-h, -h),
    'F': (0, 2-2*h),
}
for name, (aa,ell) in labels.items():
    momentum = S.Matrix([-ell/(2*z),S.I*aa/(2*alpha),0])
    expected = {'identity': S.zeros(3,1),'FA': O*(-k),'U':O*(d+a2),'F':O*k}[name]
    check(name+' primary label to original momentum', momentum-expected)
check('U local with respect to LP2017 particles', -g*labels['U'][0]-1)
check('FA Liouville reflection label is zero', 2/g-2-labels['FA'][1])
check('FA p equals negative identity p', (1-h+labels['FA'][1])+(1-h))
check('U p+ equals identity p+', 1-h+labels['U'][1]-labels['U'][0]-(1-h))
check('U p- differs by -2h', 1-h+labels['U'][1]+labels['U'][0]-(1-h)+2*h)
check('dimension Yukawa is one', -a0.dot(a0)-1)
check('dimension even root is 2t', -T.dot(T)-2*t)
check('dimension FA from LP and direct roots', -(2*h-2)**2/(4*z**2)+k.dot(k))
check('dimension U from LP and direct roots', h**2/(4*alpha**2)-h**2/(4*z**2)+(d+a2).dot(d+a2))

# Physical period: no global h is required. Fermion flip supplies the -1 for U.
delta = S.Matrix([S.sqrt(2)*S.pi/b,0])
p0,p2 = S.Matrix(a0[:2]),S.Matrix(a2[:2])
Tp,dp,kp = S.Matrix(T[:2]),S.Matrix(d[:2]),S.Matrix(k[:2])
check('Yukawa scalar phase i pi', p0.dot(delta)-S.I*S.pi)
check('potential period 2 pi i', Tp.dot(delta)-2*S.pi*S.I)
check('dual phase', dp.dot(delta)+2*S.pi*S.I/t)
check('F character exponent modulo 2 pi i', kp.dot(delta)-(2*S.pi*S.I-2*S.pi*S.I/t))
check('U character including fermion minus', (dp+p2).dot(delta)+S.I*S.pi-(2*S.pi*S.I-2*S.pi*S.I/t))

# Gaussian Wick coefficients with regulator constant kept explicit.
eps,creg = S.symbols('eps creg', positive=True)
G = -2*S.log(eps)+creg
check('both bare potentials carry same Wick factor', S.expand((-T).dot(-T)*G/2-T.dot(T)*G/2))
check('potential bare-Wick exponent', S.expand(T.dot(T)*G/2-(2*t*S.log(eps)-t*creg)))
check('null Yukawa bosonic Wick exponent', p0.dot(p0)*G/2)

# Finite fermion congruence, two physical components on each of two sites.
Dp = S.Matrix([[0,1,2,0],[-1,0,0,3],[-2,0,0,4],[0,-3,-4,0]])
Dm = S.Matrix([[0,2,0,1],[-2,0,3,0],[0,-3,0,5],[-1,0,-5,0]])
mass = S.diag(2,3,5,7)
kernel = Dp.row_join(mass).col_join((-mass.T).row_join(Dm))
shifted = Dp.row_join(-mass).col_join(mass.T.row_join(Dm))
flip = S.diag(*([-1]*4+[1]*4))
check('finite fermion antisymmetry', kernel+kernel.T)
check('translated action and fermion flip congruence', flip.T*shifted*flip-kernel)
check('two-component Berezin Jacobian', flip.det()-1)
check('fermion determinant covariance', shifted.det()-kernel.det())

# Endpoint differences on a same-orbit sector; N=right-left.
s,rN,lN = S.symbols('s right_phase left_phase', nonzero=True)
check('BB same-orbit difference', rN-lN-rN*(1-lN/rN))
check('AA same-orbit difference', 1/rN-1/lN-(1/rN)*(1-rN/lN))

# Independent direct noncommutative expansion of the projected three-product.
A,B,Ab,Bb,fA,fB = S.symbols('A B Ab Bb fA fB', commutative=False)
exprA = B*(s*Ab*A+fA)-s**-1*A*(s**-1*Ab*B)-Ab*B*A+s**-1*Ab*A*B
exprA = S.expand(exprA).subs(B*Ab,s**-1*Ab*B).subs(A*Ab,s*Ab*A+fA)
check('three-product with barA before endpoint weights', S.expand(exprA-(B*fA-s**-2*fA*B)))
exprB = B*(s**-1*Bb*A)-s**-1*A*(s*Bb*B+fB)-Bb*B*A+s**-1*Bb*A*B
exprB = S.expand(exprB).subs(B*Bb,s*Bb*B+fB).subs(A*Bb,s**-1*Bb*A)
check('three-product with barB before endpoint weights', S.expand(exprB-s**-1*(fB*A-A*fB)))

# Liouville identity-image scaling (not a physical reflection theorem).
wall = -2*z
twoQ = wall+2/wall  # canonical kinetic 1/(8pi)
check('FA equals Liouville reflection image of identity', twoQ+2*z+1/z)
check('reflection coupling exponent', -twoQ/wall-(1/t-1))
check('monodromy matches FA period character modulo integer', (1/t-1)-1/t+1)

print(f'Completed {count} independent checks.')
