"""Independent classical root construction, mass test and SG subsector test.

This tests the standard root-system super-Toda construction, NOT the different
one-boson/Dirac-fermion Toda model of Ridout--Teschner.
"""
import sympy as s
import numpy as np

def main():
    rt=s.sqrt(2); a0=s.Matrix([-1/rt,s.I/rt]);a1=s.Matrix([rt,0]);a2=s.Matrix([-1/rt,-s.I/rt])
    roots=[a0,a1,a2]
    gram=s.Matrix(3,3,lambda i,j:(roots[i].T*roots[j])[0])
    assert gram==s.Matrix([[0,-1,1],[-1,2,-1],[1,-1,0]])
    assert a0+a1+a2==s.zeros(2,1)
    print('PASS independent affine root Gram matrix:',gram.tolist())
    u,v=s.symbols('u v', real=True)
    mu0,mu1,mu2=s.symbols('mu0 mu1 mu2',positive=True)
    phi=s.Matrix([u,v]); params=[mu0,mu1,mu2]
    V=mu1*s.exp((a1.T*phi)[0])
    for i in (0,2):
        for j in (0,2):
            V+=s.Rational(1,2)*params[i]*params[j]*(roots[i].T*roots[j])[0]*s.exp(((roots[i]+roots[j]).T*phi)[0])
    V=s.simplify(V)
    assert s.simplify(V-mu1*s.exp(rt*u)-mu0*mu2*s.exp(-rt*u))==0
    H=s.hessian(V,(u,v));vac=s.log(mu0*mu2/mu1)/(2*rt)
    Hv=s.simplify(H.subs(u,vac))
    assert Hv==s.diag(4*s.sqrt(mu0*mu1*mu2),0)
    print('PASS bosonic potential=',V,'; vacuum Hessian=',Hv.tolist())
    M=s.zeros(2)
    for i in (0,2):
        M+=params[i]*roots[i]*roots[i].T*s.exp((roots[i].T*phi)[0])
    assert s.simplify(M*M-mu0*mu2*s.exp(-rt*u)*s.eye(2))==s.zeros(2)
    print('PASS fermion mass matrix squared=',s.simplify(M*M),'(conventions: normalized root fields)')
    # Independent SG kink: chi=4/k atan exp(m*x), V=m^2/k^2*(1-cos(k*chi)).
    y=s.symbols('y',positive=True)
    angle=4*s.atan(y)
    first=y*s.diff(angle,y);second=y*s.diff(first,y)
    assert s.simplify(s.trigsimp(second-s.sin(angle)))==0
    density=s.simplify(first**2/2+1-s.cos(angle))
    integral=s.integrate(s.simplify(density/y),(y,0,s.oo))
    assert integral==8
    print('PASS SG embedded kink equation; classical mass=8*m/k^2')
    # Two-kink tau solution; no use of the book scripts.
    theta=np.array([.57,-.81]);mom=np.cosh(theta);en=-np.sinh(theta)
    A=np.tanh((theta[0]-theta[1])/2)**2
    worst=0.
    for xx,tt in ((.2,.3),(-.7,.1),(1.,-.5)):
        ys=np.exp(mom*xx+en*tt)
        tau=1+1j*sum(ys)-A*np.prod(ys)
        tc=np.conj(tau)
        derivatives=[]
        for rates in (mom,en):
            d=1j*np.dot(rates,ys)-A*np.sum(rates)*np.prod(ys)
            dd=1j*np.dot(rates**2,ys)-A*np.sum(rates)**2*np.prod(ys)
            derivatives.append(-2j*((dd/tau-(d/tau)**2)-np.conj(dd/tau-(d/tau)**2)))
        Z=tau/tc;sin=(Z**2-Z**-2)/(2j)
        err=abs(derivatives[0]-derivatives[1]-sin);worst=max(worst,err)
    assert worst<1e-12,worst
    print('PASS independent SG two-kink equation; max residual',worst)
    print('SG tau interaction A=tanh^2((theta1-theta2)/2)=',A)
    print('Fast kink spatial displacement -log(A)/(m*cosh(theta1))=',-np.log(A)/mom[0])
    assert abs(np.log(A))>.1
    print('PASS same-charge SG scattering has a nonzero classical displacement; a constant S_11=1 cannot match it')

if __name__=='__main__':main()
