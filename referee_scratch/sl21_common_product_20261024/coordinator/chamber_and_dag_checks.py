"""Fresh chamber residue and whole-boundary diagnostic."""
import math
import numpy as np
from scipy.special import roots_jacobi
import sympy as s

rho=-math.sqrt(2)/5

def angular(g,n):
    # Factor endpoint powers into Jacobi weights. u^(g-2rho)(1-u)^g
    # and y^(2-rho)(1-y)^(-rho) are integrated exactly as weights.
    xu,wu=roots_jacobi(n,g,g-2*rho)
    xy,wy=roots_jacobi(n,-rho,2-rho)
    u=(xu+1)/2; y=(xy+1)/2
    wu=wu/2**(2*g-2*rho+1)
    wy=wy/2**(3-2*rho)
    return np.sum(wu[:,None]*wy[None,:]*(1+u[:,None])*(1-u[:,None]*y[None,:])**(-rho))

values=[angular(rho,n) for n in [32,64,128,256]]
for n,v in zip([32,64,128,256],values):print('INFO chamber angular',n,format(v,'.14g'),'residue/3',format(v/3,'.14g'))
assert values[-1]>0 and abs(values[-1]-values[-2])<2e-8
print('PASS nonzero individual-chamber residue independently integrated')
for delta in [0.02,0.005,0.001]:
    # delta times Mellin radial integral is exactly 1/3 for Re delta>0.
    print('INFO delta',delta,'finite residue',angular(rho+delta,128)/3)

R,u,y=s.symbols('R u y',positive=True)
coords=s.Matrix([R,u*R,u*y*R])
det=coords.jacobian([R,u,y]).det()
assert s.factor(det-u*R**2)==0
print('PASS full three-relative-coordinate Jacobian u R^2')

# DAG properness proof in finite dimensions: for oriented edge i->j,
# i<j, choose H_i=n-i. Then H.(incidence k)=sum_(i<j)(j-i)kij.
# Each coefficient >=1 bounds total Schwinger parameters on a fixed image.
for n in range(2,7):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    H=np.arange(n,0,-1)
    for i,j in edges:
        e=np.zeros(n,dtype=int);e[i]=1;e[j]=-1
        assert np.dot(H,e)==j-i>=1
print('PASS acyclic incidence has a strictly positive linear functional')
print('INFO Thus positive edge parameters cannot run to infinity at fixed momentum; normalized power distributions push forward without a positive cycle.')

# A two-edge opposite orientation has a nonzero positive circulation.
inc=np.array([[1,-1],[-1,1]])
assert np.array_equal(inc@np.array([1,1]),[0,0])
print('PASS holomorphic-antiholomorphic pair gives a positive cycle, outside the DAG proof')

# The simplest exact demonstration that chamber residues need not survive:
# boundary power x_+^lambda + exp(i*pi*lambda)(-x)_+^lambda has cancelling
# first Taylor residues at lambda=-2, unlike either chamber separately.
L=s.symbols('L')
assert s.simplify((1-s.exp(s.I*s.pi*L)).subs(L,-2))==0
assert s.simplify((1+s.exp(s.I*s.pi*L)).subs(L,-1))==0
print('PASS both first two one-dimensional chamber residues cancel in complete i0 boundary value')
print('INFO No full five-current anomaly is inferred from the nonzero chamber number.')
