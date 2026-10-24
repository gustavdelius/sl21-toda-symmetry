"""Independent dual-screening OPE calculation, canonical UV free bosons.
X X ~ -log; coefficients in exponent are complex momenta, no extra i.
Fermions chi0=exp(i h), chi2=exp(-i h), with suitable Klein factors.
All claims are first-order unless explicitly marked as an obstruction.
"""
import sympy as s
b=s.symbols('b',nonzero=True)
r=s.sqrt(2)
a0=s.Matrix([s.I*b/r,b/r,s.I]);a2=s.Matrix([s.I*b/r,-b/r,-s.I]);a1=-a0-a2
T=a0+a2
D=lambda x,y:s.expand((x.T*y)[0])
assert D(a0,a0)==D(a2,a2)==-1
assert D(a0,a2)==1-b*b
assert D(T,T)==-2*b*b
print('PASS bosonization gives fermionic screening norms -1 and pair product 1-b^2')
d=(2*T/D(T,T)).applyfunc(s.simplify);c=-d
assert d==s.Matrix([-s.I*r/b,0,0])
assert D(d,a0)==D(d,a2)==1 and D(d,a1)==-2
print('PASS pair-dual momentum is -i*alpha1/b, opposite to the even-node BL current')
# For J=(C.dX)e^(dX) and V=e^(aX), d.a=1, residue is
# [C+(C.a)*a].dX e^((d+a)X).
C=a0
R0=C+D(C,a0)*a0
R2=C+D(C,a2)*a2
assert R0==s.zeros(3,1)
assert (R2+b*b*(d+a2)).applyfunc(s.simplify)==s.zeros(3,1)
print('PASS node-0 residue vanishes; node-2 residue = -b^2*d_z exp((d+a2).X)')
print('PASS even-root perturbation has regular OPE with the dressed dual current')
# Even-node J+ has double pole with a1, regular OPEs with a0,a2,T.
gamma=b*b/(1-b*b)
assert (a1+gamma*(c+a1)).applyfunc(s.simplify)==s.zeros(3,1)
assert D(c,a0)==D(c,a2)==-1 and D(c,T)==-2
print('PASS J+ residue reproduces -gamma*d_z exp((c+a1).X), other first-order OPEs regular')
# Potential quantum-group braiding data, not Serre or charge-algebra proof.
M=s.Matrix([[D(c,c),D(c,d)],[D(d,c),D(d,d)]])
assert (M/(-1/b**2))==s.Matrix([[2,-2],[-2,2]])
print('PASS exponential momentum Gram matrix is (-1/b^2)*A1^(1) Cartan matrix; both currents are even')
# Extra free direction, not used to claim exact decoupling of full quantum model.
R=s.Matrix([0,r,s.I*b])
assert D(R,a0)==D(R,a2)==D(R,a1)==0 and D(R,R)==2-b*b
print('PASS an orthogonal free-field direction exists in the formal bosonized perturbation system')
# Composite bosonic perturbation V_T: d.T=2. Its residue includes level-two
# jets. Subtract an explicit total derivative and exhibit remainder.
x=s.Matrix(s.symbols('x1:4'));xx=s.Matrix(s.symbols('xx1:4'))
CA=D(C,T);P=d+T
res=D(C,x)*D(T,x)+CA*(D(T,xx)+D(T,x)**2)/2
hder=CA*(D(T,xx)+D(T,x)*D(P,x))/2
remainder=D(T,x)*D(a0-a2,x)/2
assert s.simplify(res-hder-remainder)==0
assert s.simplify(remainder)!=0
print('PASS composite-potential residue modulo derivatives has explicit nonzero representative 1/2*(T.dX)*((a0-a2).dX)*exp((d+T).X)')
print('CAUTION a double-Yukawa contribution occurs at the same mass order; cancellation is NOT established')
# Spins with zero background charge: vertex h=-momentum^2/2.
hplus=-D(c,c)/2;hminus=1-D(d,d)/2
assert s.simplify(hplus-1-(1/b**2-1))==0
assert s.simplify(hminus-1-1/b**2)==0
print('PASS free-UV charge spins are 1/b^2-1 and 1/b^2; not fixed physical spins before renormalization')
