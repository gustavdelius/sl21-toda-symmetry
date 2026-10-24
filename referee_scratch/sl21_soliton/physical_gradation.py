"""Make categorical crossing into an explicit rapidity crossing convention.

G_v(theta)=diag(exp(-k theta),1,exp(theta/2)), G_b=G_v^-1.
This is a chosen gradation, not derived spins of Toda nonlocal currents.
"""
from candidate import *

def st1n(a):
    b=np.zeros((9,9),complex)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    b[3*k+j,3*i+l]=(-1)**((p[j]+p[l])*(p[k]+p[i])+p[i]*(p[i]+p[k]))*a[3*i+j,3*k+l]
    return b

def main():
    mats=physical_candidate()
    f={n:s.lambdify((t,x),r,'numpy') for n,r in mats.items()}
    worst=0.
    for kap in (np.sqrt(2)/3,np.sqrt(3)/2, np.sqrt(2)):
        tn=np.exp(.5j*np.pi*kap);qn=tn**2
        C=np.diag([qn**-1,1,1j]);C1=np.kron(C,np.eye(3))
        def g(a,theta):
            v=np.array([np.exp(-kap*theta),1,np.exp(theta/2)],complex)
            return np.diag(v if a=='v' else 1/v)
        def r(a,b,theta):
            G=np.kron(g(a,theta),np.eye(3))
            return G@f[a+b](tn,np.exp(kap*theta))@np.linalg.inv(G)
        for a in 'vb':
            for b in 'vb':
                theta=.39+.07j
                crossed=np.linalg.inv(C1)@st1n(np.linalg.inv(r(a,b,theta)))@C1
                target=r('b' if a=='v' else 'v',b,theta+1j*np.pi)
                err=np.max(np.abs(crossed-target));worst=max(worst,err)
                assert err<2e-12,(a,b,kap,err)
        # Unitarity uses the actual two-rapidity matrices, not a guessed fixed gauge.
        for a in 'vb':
            for b in 'vb':
                # R_ba(-theta), conjugated through swap, has the equivalent first
                # rapidity gauge; total Cartan covariance supplies the identity.
                err=np.max(np.abs(r(a,b,.39)@np.array(P,float)@r(b,a,-.39)@np.array(P,float)-np.eye(9)))
                assert err<1e-11,err
        sv=np.linalg.svd(np.array(P,float)@r('v','v',.39),compute_uv=False)
        print('kappa',kap,'crossing-gradation vv singular values',sv)
    print('PASS all 4 crossed channels and braiding unitarity, 3 generic kappas; crossing max residual',worst)
    print('Chosen generator spins (e0,e1,e2)=(1/2,-kappa,-1/2), not a coupling dictionary.')

if __name__=='__main__':main()
