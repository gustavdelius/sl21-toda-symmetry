"""Structural strip-pole enumeration and independent 1/epsilon confirmation."""
from candidate import *

def locations(kappa, phase):
    return [(phase+2*np.pi*n)/kappa for n in range(-20,21)
            if 1e-10<(phase+2*np.pi*n)/kappa<np.pi-1e-10]

def main():
    phys=physical_candidate()
    f={n:s.lambdify((t,x),r,'numpy') for n,r in phys.items()}
    for kappa in (np.sqrt(2)/3,np.sqrt(3)/2,np.sqrt(2),np.sqrt(3)):
        tn=np.exp(.5j*np.pi*kappa)
        vv=locations(kappa,2*np.pi*kappa)
        mixed=locations(kappa,-np.pi*kappa)
        print('kappa=',kappa,'VV poles u/pi=',[a/np.pi for a in vv],
              'mixed poles u/pi=',[a/np.pi for a in mixed])
        if 0<kappa<1: assert not vv and not mixed
        if 1<kappa<2:
            assert len(vv)==len(mixed)==1
            assert abs(vv[0]+mixed[0]-np.pi)<1e-12
            for name,us in (('vv',vv),('bv',mixed),('vb',mixed)):
                theta=1j*us[0]
                norms=[np.linalg.norm(f[name](tn,np.exp(kappa*(theta+eps)))) for eps in (1e-5,1e-6)]
                ratio=norms[1]/norms[0]
                assert abs(ratio-10)<1e-3,(name,ratio)
                print('PASS structural simple pole confirmed by 1/epsilon scaling',name,'ratio=',ratio)
    # Explicit scalar ambiguity: this factor is unitary and crossing-symmetric,
    # independently of the matrix part. a is an arbitrary CDD parameter.
    theta,a=s.symbols('theta a',real=True)
    h=(s.sinh(theta)+s.I*s.sin(a))/(s.sinh(theta)-s.I*s.sin(a))
    assert s.simplify(h*h.subs(theta,-theta))==1
    assert s.simplify(h.subs(theta,s.I*s.pi-theta)-h)==0
    assert s.simplify(h.subs(theta,theta+2*s.I*s.pi)-h)==0
    print('PASS scalar CDD freedom: unitarity, crossing and 2pi*i periodicity')

if __name__=='__main__':main()
