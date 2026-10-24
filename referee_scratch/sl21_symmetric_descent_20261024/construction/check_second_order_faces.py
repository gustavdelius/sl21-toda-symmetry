import sympy as s
from itertools import product
n=0
for equal in [False,True]:
 for a,b,d in product(range(2),range(4),range(2)):
  if a+b+d==0:continue
  E=7*b*(b-1)/2+2*a*b
  degree=E+2*a-2*b-2*d+(5 if equal else 2*8+1)
  assert degree>0,(equal,a,b,d,degree)
  print('PASS '+('equal' if equal else 'unlike')+' two-source union '+str((a,b,d))+' degree='+str(degree));n+=1
# Formal direct expansion of once-per-original-slot replacement.
l0,l2=s.symbols('l0 l2')
f=s.symbols('f0:5'); u0=s.symbols('u00:5');u2=s.symbols('u20:5')
p=s.prod(f[j]+l0*u0[j]+l2*u2[j] for j in range(5))
direct=s.expand(p).coeff(l0,1).coeff(l2,1)
forest=sum(u0[i]*u2[j]*s.prod(f[k] for k in range(5) if k not in (i,j)) for i in range(5) for j in range(5) if i!=j)
assert s.expand(direct-forest)==0
print('PASS exact distinct original-slot U02 combinatorics');n+=1
for i in range(5):
 assert not any(s.degree(term,u0[i]) and s.degree(term,u2[i]) for term in s.Add.make_args(direct))
print('PASS no same-original-slot double replacement');n+=1
# Retain order-labelled spatial insertions as independent formal linear maps.
S0,S2,S02,F,U0,U2,U02=s.symbols('S0 S2 S02 F U0 U2 U02',commutative=False)
# Merely checks coefficient weights; it does not identify S0 U2 with U2 S0.
P=(1+l0*S0+l2*S2+l0*l2*S02)*(F+l0*U0+l2*U2+l0*l2*U02)
coeff=s.expand(P).coeff(l0,1).coeff(l2,1)
assert s.expand(coeff-(S02*F+S0*U2+S2*U0+U02))==0
print('PASS insertion/repair order02 coefficient preserves map order');n+=1
print('PASS total '+str(n))
