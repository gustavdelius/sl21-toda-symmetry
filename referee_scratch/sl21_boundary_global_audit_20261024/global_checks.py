"""Fresh local/canonical/finite-cutoff audit; imports no previous check code."""
import sympy as sp
import mpmath as mp

passed = 0
def check(label, value):
    global passed
    if value is not True:
        raise AssertionError('FAIL ' + label + ': ' + str(value))
    passed += 1
    print('PASS ' + label)

def zero(label, expr):
    check(label, sp.simplify(expr) == 0)

b = sp.symbols('b', nonzero=True, real=True)
t = b*b
I = sp.I
p0 = sp.Matrix([I*b/sp.sqrt(2), b/sp.sqrt(2)])
p2 = sp.Matrix([I*b/sp.sqrt(2), -b/sp.sqrt(2)])
T = p0+p2
d = -T/t
k = d+T
delta = sp.Matrix([sp.sqrt(2)*sp.pi/b, 0])
for name, p in [('0',p0), ('2',p2)]:
    zero('physical null root '+name, p.dot(p))
    zero('dual pairing '+name, d.dot(p)-1)
    zero('Yukawa scalar period '+name, p.dot(delta)-I*sp.pi)
zero('T period', T.dot(delta)-2*I*sp.pi)
zero('dual period', d.dot(delta)+2*I*sp.pi/t)
zero('F phase differs from s by 2 pi i', k.dot(delta)+2*I*sp.pi/t-2*I*sp.pi)
zero('U phase includes direct fermion minus', (d+p2).dot(delta)+I*sp.pi+2*I*sp.pi/t-2*I*sp.pi)
zero('FA inverse phase', -k.dot(delta)-2*I*sp.pi/t+2*I*sp.pi)
zero('dual length', d.dot(d)+2/t)
zero('topological dual pairing', T.dot(d)-2)

# Bare lattice fermionic kinetic blocks are invariant under two left sign flips.
D = sp.diag(-1,-1,1,1)
m0,m2 = sp.symbols('m0 m2')
A = sp.Matrix([[0,0,I*m0,0],[0,0,0,I*m2],[-I*m0,0,0,0],[0,-I*m2,0,0]])
def pf4(M): return M[0,1]*M[2,3]-M[0,2]*M[1,3]+M[0,3]*M[1,2]
zero('two-component left Jacobian', D.det()-1)
check('mass sign change from left parity', D.T*A*D == -A)
zero('single-site mass Pfaffian', pf4(A)-m0*m2)
zero('Pfaffian unchanged under left parity', pf4(D.T*A*D)-pf4(A))
check('left kinetic bilinears unchanged', sp.diag(-1,-1).T*sp.eye(2)*sp.diag(-1,-1)==sp.eye(2))

# General finite anchor convention. eta=0 is the report's left anchor;
# eta=1 is a right anchor. Only N's commutator is anchor independent.
eta,q = sp.symbols('eta q')
left_weight = (eta-1)*q*T.dot(d)
right_weight = eta*q*T.dot(d)
zero('topological weight independent of anchor', right_weight-left_weight-2*q)
zero('left anchored right endpoint fixed', right_weight.subs(eta,0))
zero('left anchored left endpoint shift', left_weight.subs(eta,0)+2*q)
zero('right anchored left endpoint fixed', left_weight.subs(eta,1))
zero('right anchored right endpoint shift', right_weight.subs(eta,1)-2*q)
tails = [('A',-d,sp.zeros(2,1),-1),('B',d,sp.zeros(2,1),1),
         ('HA',-k,-T,-1),('HB1',d+p2,p2,1),('HB2',k,T,1)]
for name,l,r,sign in tails:
    check(name+' has same disorder momentum as its current', sp.simplify(l-r-sign*d)==sp.zeros(2,1))
    zero(name+' topological weight', T.dot(l-r)-2*sign)
zero('F left BCH phase versus character', -2*I*sp.pi*k.dot(d)-4*I*sp.pi/t+4*I*sp.pi)
zero('U left BCH phase versus character', -2*I*sp.pi*(d+p2).dot(d)-4*I*sp.pi/t+2*I*sp.pi)

# EOM projection coefficients are computed directly from physical roots.
for name,p,expected in [('Y0',p0,1),('Y2',p2,1),('V1',-T,-2),('P',T,2)]:
    zero('dual EOM coefficient '+name,d.dot(p)-expected)
S=p0-p2
zero('difference EOM selects Y0',S.dot(p0)-t)
zero('difference EOM selects minus Y2',S.dot(p2)+t)
zero('difference EOM removes both potentials',S.dot(T))

# Changing the Gaussian subtraction constant rescales both bosonic potentials.
eps,c = sp.symbols('eps c', positive=True)
G=-2*sp.log(eps)+c
for name,p in [('P',T),('V1',-T)]:
    zero(name+' bare/Wick logarithmic conversion',p.dot(p)*G/2-(2*t*sp.log(eps)-t*c))
for name,p in [('Y0 scalar',p0),('Y2 scalar',p2)]:
    zero(name+' has no bosonic self factor',p.dot(p)*G/2)

# Exact one-dimensional integration-by-parts probe. This is not Toda QFT.
# It checks the contact term in the weighted finite-cutoff SD identity.
mp.mp.dps=50
a=mp.sqrt(2); bb=mp.e/2; cnum=mp.sqrt(3)/2
def moment(power):
    return 2*(bb/a)**(power/2)*mp.besselk(power,2*mp.sqrt(a*bb))
res=a*moment(cnum+1)-bb*moment(cnum-1)-cnum*moment(cnum)
check('weighted SD retains derivative of insertion',bool(abs(res)<mp.mpf('1e-45')))
res_wrong=a*moment(cnum+1)-bb*moment(cnum-1)
check('discarding weighted SD contact is false',bool(abs(res_wrong)>mp.mpf('1e-4')))
quad=mp.quad(lambda x:mp.exp(cnum*x-a*mp.exp(x)-bb*mp.exp(-x)),[-12,0,12])
check('SD probe Bessel moment independently integrated',bool(abs(quad-moment(cnum))<mp.mpf('1e-40')))

print(f'TOTAL {passed} PASS, 0 FAIL')
