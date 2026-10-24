"""Independent UV OPE tests in Bernard-LeClair canonical normalization.
phi^a phi^b ~ -delta^{ab} log(w-z), psi^a psi^b ~ delta^{ab}/(w-z).
This is first-order conformal perturbation theory, not an all-order claim.
"""
import sympy as s
b=s.symbols('b',nonzero=True)
r=s.sqrt(2)
roots=[s.Matrix([-1/r,s.I/r]),s.Matrix([r,0]),s.Matrix([-1/r,-s.I/r])]
G=s.Matrix([[(a.T*c)[0] for c in roots] for a in roots])
assert G==s.Matrix([[0,-1,1],[-1,2,-1],[1,-1,0]])
# V_1(w)J_1(z) has double pole and residue (-i b alpha1.dphi) E.
# E=exp(i(1/b-b) alpha1.phi). Gamma=b^2/(1-b^2).
alpha=roots[1];gamma=b*b/(1-b*b)
assert (-(gamma)*s.I*(1/b-b)*alpha+s.I*b*alpha).applyfunc(s.simplify)==s.zeros(2,1)
print('PASS even-root BL residue is -gamma*d_z exp(i(1/b-b)*alpha1.phi), gamma=b^2/(1-b^2)')
assert G[1,0]==G[1,2]==-1
assert (roots[1].T*(roots[0]+roots[2]))[0]==-2
print('PASS even-root vertex has regular OPE with both Yukawa holomorphic vertices and the composite bosonic exponential')
# Bare odd ansatz J=(p.psi)exp(i alpha0.phi/b).
# Yukawa node 0: simple pole (alpha0.p) E, requiring alpha0.p=0.
# Yukawa node 2: residue contains :(alpha2.psi)(p.psi): E;
# its coefficient is det(alpha2,p), absent from a derivative of a plain vertex.
p0,p1=s.symbols('p0 p1');p=s.Matrix([p0,p1])
no_simple=(roots[0].T*p)[0]
no_bilinear=s.det(s.Matrix.hstack(roots[2],p))
sol=s.solve([no_simple,no_bilinear],(p0,p1),dict=True)
assert sol==[{p0:0,p1:0}]
print('PASS no nonzero bare linear-fermion prefactor solves both odd-node tests with plain-vertex counterterms')
print('CAUTION derivative prefactors, additional disorder dressing and other UV realizations remain untested; no no-go theorem follows')
