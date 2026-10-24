import sympy as s
l0,l2=s.symbols('l0 l2')
f=s.symbols('f0:5'); u0=s.symbols('u0_0:5');u2=s.symbols('u2_0:5')
p=s.prod(f[j]+l0*u0[j]+l2*u2[j] for j in range(5))
direct=s.expand(p).coeff(l0,1).coeff(l2,1)
forest=sum(u0[i]*u2[j]*s.prod(f[k] for k in range(5) if k not in (i,j)) for i in range(5) for j in range(5) if i!=j)
assert s.expand(direct-forest)==0
print('PASS exact distinct original-slot U02 combinatorics')
for i in range(5):
 assert not any(s.degree(term,u0[i]) and s.degree(term,u2[i]) for term in s.Add.make_args(direct))
print('PASS no same-original-slot double replacement')
S0,S2,S02,F,U0,U2,U02=s.symbols('S0 S2 S02 F U0 U2 U02',commutative=False)
P=(1+l0*S0+l2*S2+l0*l2*S02)*(F+l0*U0+l2*U2+l0*l2*U02)
coeff=s.expand(P).coeff(l0,1).coeff(l2,1)
assert s.expand(coeff-(S02*F+S0*U2+S2*U0+U02))==0
print('PASS insertion/repair order02 coefficient preserves map order')
print('PASS total 3')
