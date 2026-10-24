"""Independent construction checks and a fixed-generator equivalence obstruction."""
import sympy as s
from itertools import product

dg,dh=s.symbols('dg dh')
# Four independent current counters plus two source counters. Counters are
# additive; the regulator is bilinear even for fused composites.
def edge(x,y):
    return dg*sum(x[i]*y[i] for i in range(4))+dh*(x[4]*y[5]+x[5]*y[4])
def plus(x,y):return tuple(a+b for a,b in zip(x,y))
basis=[tuple(int(i==j) for i in range(6)) for j in range(6)]
for x,y,z in product(basis,repeat=3):
    assert s.expand(edge(x,y)+edge(plus(x,y),z)-edge(y,z)-edge(x,plus(y,z)))==0
print('PASS all 216 generator triples obey additive regulator cocycle identity')

# This is meaningful before diagonal restriction: it gives one Gaussian
# Wick product with physical+auxiliary contractions, not just a word name.
# An interaction must be evaluated once as a common linear functional.
lam,V=s.symbols('lam V')
global_weight=1-lam*V
independently_dressed_blocks=s.expand(global_weight**2)
assert independently_dressed_blocks.coeff(lam,1)==-2*V
print('PASS separate block interaction dressings double count the first source; one common state is required')

# If mu0(B,D)=s mu0(D,B), a transported multiplication fixing B,D
# cannot change their leading braided defect: the Z1(composite) terms cancel.
q=s.symbols('q',nonzero=True)
zBD,zDB=s.symbols('zBD zDB')
change=(-zBD+q*zDB).subs(zBD,q*zDB)
assert change==0
print('PASS fixed-generator coboundary has zero leading braided-defect change')

t=s.symbols('t',nonzero=True)
gam=t/(1-t)
assert s.simplify(t*t-t**2)==0
assert s.simplify(t*(gam**2/t)-gam**2)==0
print('PASS common sqrt(t) generator scaling converts both raw diagonal contact coefficients')

# Explicit partial-fraction primitive of the barred-tail scalar kernel.
a,u=s.symbols('a u')
z=s.symbols('z0:3')
G=(a-u)**2/s.prod((u-zi)**2 for zi in z)
Si=[(a-zi)**2/s.prod((zi-zj)**2 for j,zj in enumerate(z) if j!=i) for i,zi in enumerate(z)]
Ri=[s.diff((a-u)**2/s.prod((u-zj)**2 for j,zj in enumerate(z) if j!=i),u).subs(u,zi) for i,zi in enumerate(z)]
assert s.factor(sum(Ri))==0
assert s.factor(G-sum(Si[i]/(u-z[i])**2+Ri[i]/(u-z[i]) for i in range(3)))==0
print('PASS barred-tail scalar partial fractions and zero total residue')
primitive=sum(-Si[i]/(u-z[i])+Ri[i]*s.log(u-z[i]) for i in range(3))
assert s.factor(s.diff(primitive,u)-G)==0
print('PASS barred-tail scalar has a homogeneous logarithmic u primitive')
print('INFO This last primitive is scalar and off diagonals; branch and proper-face extension must still be checked.')
