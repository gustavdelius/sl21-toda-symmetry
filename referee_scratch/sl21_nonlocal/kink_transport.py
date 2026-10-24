"""First Grassmann derivative of vacuum-subtracted transport on the exact kink.
This tests convergence, not charge conservation; finite cutoffs have nonzero flux.
"""
import numpy as np
from scipy.integrate import quad_vec
from scipy.linalg import expm
mu=np.sqrt(2)
B=np.array([[0,1j*mu,0],[-1j*mu,0,0],[0,0,0]],complex)
root=[np.array([-1,1j])/np.sqrt(2),np.array([-1,-1j])/np.sqrt(2)]
def insertion(x,which):
 y=np.exp(2*mu*x);f=np.sqrt(y)/(1+1j*y)
 ss=(1-1j*y)/np.sqrt(1+y*y)
 pp=np.array([f if which==0 else 0, f if which==1 else 0]);pm=pp*np.array([1,-1])
 a0p,a2p=[r@pp for r in root];a0m,a2m=[r@pm for r in root]
 mat=np.zeros((3,3),complex)
 mat[2,0]=-1j*np.sqrt(mu)*ss*a0p/2
 mat[1,2]=-1j*np.sqrt(mu)*ss*a2p/2
 mat[0,2]=np.sqrt(mu)*ss*a0m/2
 mat[2,1]=-np.sqrt(mu)*ss*a2m/2
 return mat
# Body matrices are even, so no extra Koszul sign enters the first variation.
def integrand(x,k):
 h=np.log(np.cosh(2*mu*x))/(2*mu)
 return -expm(-h*B)@insertion(x,k)@expm(h*B)
vals=[]
for ell in (2.,4.,6.,8.):
 L=ell/mu;h=np.log(np.cosh(2*mu*L))/(2*mu)
 left=expm((-L+h)*B);right=expm((L-h)*B)
 norms=[]
 for k in (0,1):
  integral,err=quad_vec(lambda x:integrand(x,k),-L,L,epsabs=1e-10,epsrel=1e-10)
  result=left@integral@right
  norms.append(np.linalg.norm(result))
 vals.append(norms)
 print('CUTOFF mu*L',ell,'odd first-variation norms',norms,flush=True)
expected=1/(2*np.sqrt(2*mu))
assert abs(vals[-1][0]-expected)<1e-9
assert abs(vals[-1][1]-expected)<1e-9
print('PASS both odd first variations converge to the independently derived norm 1/(2 sqrt(2 mu)) =',expected)
