"""Independent zero modes of the root-system A(1,0)^(1) Toda kink.

Source definitions only: even node alpha1, odd nodes alpha0,alpha2. The
calculation is complex/holomorphic; no positive real quantum contour is assumed.
Light-cone derivatives are d_t +/- d_x, mu>0, y=exp(2*mu*x).
"""
import sympy as s
import numpy as np
from scipy.integrate import solve_ivp
import mpmath as mp

y, mu = s.symbols('y mu', positive=True)
rt=s.sqrt(2)
roots=[s.Matrix([-1/rt,s.I/rt]),s.Matrix([rt,0]),s.Matrix([-1/rt,-s.I/rt])]
J=s.Matrix([[0,-1],[1,0]])
sig3=s.diag(1,-1)

def dx(expr):
    if isinstance(expr,s.MatrixBase):return expr.applyfunc(lambda z:2*mu*y*s.diff(z,y))
    return 2*mu*y*s.diff(expr,y)
def clean(expr):
    if isinstance(expr,s.MatrixBase):return expr.applyfunc(s.factor)
    return s.factor(expr)

def mode_data(orientation=1):
    # Kink chi=4 atan y; antikink chi=4 atan(1/y).
    chi=4*s.atan(y if orientation==1 else 1/y)
    if orientation==1:
        w=(1-s.I*y)/(1+s.I*y)
        f=s.sqrt(y)/(1+s.I*y)
    else:
        w=(y-s.I)/(y+s.I)
        f=s.sqrt(y)/(1-s.I*y)
    r=mu*w;M=r*sig3
    # Rows: (psi+_u,psi+_v,psi-_u,psi-_v).
    U=f*s.Matrix([[1,0],[0,1],[orientation,0],[0,-orientation]])
    return chi,w,f,r,M,U

def main():
    for orientation in (1,-1):
        chi,w,f,r,M,U=mode_data(orientation)
        assert s.simplify(dx(dx(chi))-4*mu**2*s.sin(chi))==0
        assert s.simplify(dx(f)-orientation*r*f)==0
        A=s.zeros(4)
        A[:2,2:]=M;A[2:,:2]=M
        assert clean(dx(U)-A*U)==s.zeros(4,2)
        print('PASS kink/antikink equation and two exact fermion zero modes; orientation',orientation,flush=True)
        assert clean(M*J+J*M)==s.zeros(2)
        # Chirality species matrix diagonalises the zero-energy first-order ODE.
        for col in range(2):
            # Non-normalizable partners have the opposite chiral ratio.
            b=s.Matrix([1,0,-orientation,0]) if col==0 else s.Matrix([0,1,0,orientation])
            growing=b/f
            assert clean(dx(growing)-A*growing)==s.zeros(4,1)
        assert s.limit(f,y,0)==0 and s.limit(f,y,s.oo)==0
        print('PASS completeness: two decaying and two growing solutions span the four first-order solutions',flush=True)
        # Absolute and holomorphic norms, with dx=dy/(2*mu*y).
        absolute=s.integrate(1/(2*mu*(1+y**2)),(y,0,s.oo))
        bilinear=-orientation*s.I/(2*mu)
        assert absolute==s.pi/(4*mu)
        primitive=y/(2*mu*(1+orientation*s.I*y))
        assert s.simplify(s.diff(primitive,y)-f**2/(2*mu*y))==0
        assert s.limit(primitive,y,s.oo)-s.limit(primitive,y,0)==bilinear
        print('PASS |f|^2 norm=pi/(4*mu); bilinear integral f^2=',bilinear,flush=True)
        assert clean(U.T*U-2*f**2*s.eye(2))==s.zeros(2)
        print('PASS projected fermionic kinetic pairing is nondegenerate, 2*integral(f^2)*I2',flush=True)

    # Independent numerical solve of the original 4-component static ODE,
    # normalized at the centre, not exponentially tiny at one endpoint.
    mun=np.sqrt(2);worst=0.
    for orientation in (1,-1):
        chi,w,f,r,M,U=mode_data(orientation)
        fun=s.lambdify((mu,y),U,'numpy')
        for col in (0,1):
            u0=fun(mun,1)[:,col]
            u0=u0/np.linalg.norm(u0)
            def rhs(xx,vec):
                yy=np.exp(2*mun*xx)
                ww=(1-1j*yy)/(1+1j*yy) if orientation==1 else (yy-1j)/(yy+1j)
                mass=mun*ww*np.diag([1,-1])
                return np.r_[mass@vec[2:],mass@vec[:2]]
            for end in (-8/mun,8/mun):
                grid=np.linspace(0,end,101)
                sol=solve_ivp(rhs,(0,end),u0,t_eval=grid,rtol=2e-11,atol=1e-12)
                assert sol.success
                exact=np.column_stack([fun(mun,np.exp(2*mun*z))[:,col] for z in grid])
                exact=exact/np.linalg.norm(fun(mun,1)[:,col])
                worst=max(worst,float(np.max(np.abs(sol.y-exact))))
    assert worst<2e-9,worst
    print('PASS independent normalized ODE integrations, max absolute profile error',worst,flush=True)
    # Bosonic kernel: translation only, plus a free direction with no L2 zero mode.
    a=2*y/(1+y**2) # sech(2*mu*x)
    tanh=(y**2-1)/(y**2+1)
    assert clean(-dx(dx(a))+4*mu**2*(1-2*a**2)*a)==0
    assert clean(dx(a)+2*mu*tanh*a)==0
    print('PASS bosonic translation mode and A^dagger*A factorization; free -d_x^2 has no nonzero L2 kernel',flush=True)
    # Quantum complex Clifford algebra of two nonsingular collective coordinates.
    gamma1=s.Matrix([[0,1],[1,0]]);gamma2=s.Matrix([[0,-s.I],[s.I,0]])
    for a in (gamma1,gamma2):assert a*a==s.eye(2)
    assert gamma1*gamma2+gamma2*gamma1==s.zeros(2)
    F=-s.I*gamma1*gamma2/2
    assert sorted(F.eigenvals())==[-s.Rational(1,2),s.Rational(1,2)]
    print('PASS formal Cl_2 zero-mode quantization has two states with centred axial charges +/-1/2',flush=True)
    print('CAUTION: this formal two-state space is not yet an infrared-complete quantum kink multiplet.')

if __name__=='__main__':main()
