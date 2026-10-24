import sympy as s
p0,pz,pw,pwb,C,Cb,beta,sg=s.symbols('psi0 psiz psiw psiwb C Cbar beta sigma')
M=s.pi*sg; J=s.I*s.pi*sg
outer=-s.I*s.pi*beta*C*p0*M+beta**2*pz*(-s.I*s.pi*M+s.pi*J)+s.pi*beta**2*(pw+pwb)*J-s.pi*beta*Cb*p0*J
expected=s.I*s.pi**2*sg*(beta**2*(pw+pwb)-beta*(C+Cb)*p0)
assert s.expand(outer-expected)==0
print('PASS complete independent mixed first-jet functional')
assert s.diff(outer,pz)==0
print('PASS arbitrary outer-slot first derivative cancels')
# Derivatives at the two unary source allocations differ only by source jets.
r=s.symbols('r'); A=p0+(pz+pw+pwb)*r; B=p0+pz*r
assert s.diff(A-B,r)==pw+pwb
print('PASS unary doublepole jump produces exactly the source jet')
# delta-prime jump acts on A-B, so zero zeroth coefficient, no outer-test jet.
assert (A-B).subs(r,0)==0
print('PASS unary delta-prime leading coefficient cancels')
naive=outer.subs({p0:0,pz:0,pw:1,pwb:0})
assert s.simplify(naive-s.I*s.pi**2*sg*beta**2)==0
print('FAIL EXPECTED constant-spectator mixed contact lift: misses i*pi^2*sigma*beta^2')
# AA/V1 has no derivative current, but the same independent source jet.
aa=-s.I*s.pi*pz*M+s.pi*(pz+pw+pwb)*J
assert s.expand(aa-s.I*s.pi**2*sg*(pw+pwb))==0
print('PASS complete independent AA/V1 first-jet functional')
print('PASS total 5; retained expected failed shortcut 1')
