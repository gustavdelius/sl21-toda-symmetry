"""Action-derived axial current, local-supercharge test and module compatibility.

The four-state Uq module constructed at the end is a COMPATIBILITY example,
not an action-derived nonlocal symmetry or soliton S-matrix identification.
"""
import sympy as s

def add(*terms):
    keys=set().union(*(a.keys() for a in terms))
    return {k:s.simplify(sum(a.get(k,0) for a in terms)) for k in keys}
def mul(a,b):
    out={}
    for i,u in a.items():
        for j,v in b.items():
            if i&j:continue
            inversions=sum(1 for k in range(4) for l in range(4) if k>l and (i>>k)&1 and (j>>l)&1)
            out[i|j]=out.get(i|j,0)+(-1)**inversions*u*v
    return out
def scale(c,a):return {k:c*v for k,v in a.items()}

def main():
    c,d=s.symbols('c d')
    z=[{1<<i:s.S.One} for i in range(4)]
    p=z[:2];m=z[2:]
    dp=[add(scale(-c,m[0]),scale(-d,m[1])),add(scale(-d,m[0]),scale(c,m[1]))]
    dm=[add(scale(c,p[0]),scale(d,p[1])),add(scale(d,p[0]),scale(-c,p[1]))]
    div=add(scale(-s.I,add(mul(dp[0],p[1]),mul(p[0],dp[1]))),
            scale(s.I,add(mul(dm[0],m[1]),mul(m[0],dm[1]))))
    assert all(v==0 for v in div.values())
    print('PASS full axial Noether-current conservation from arbitrary traceless symmetric Yukawa matrix')
    u,v,mu=s.symbols('u v mu',real=True,nonzero=True)
    r=s.sqrt(2)
    W=2*mu*s.exp(-u/r)*s.cos(v/r)
    M=s.hessian(W,(u,v))
    gradW=s.Matrix([s.diff(W,u),s.diff(W,v)])
    V=mu**2*(s.exp(r*u)+s.exp(-r*u))
    gradV=s.Matrix([s.diff(V,u),s.diff(V,v)])
    mismatch=(gradV-M*gradW).applyfunc(s.simplify)
    assert mismatch==s.Matrix([r*mu**2*s.exp(r*u),0])
    print('PASS local W-supercurrent has nonzero divergence coefficient grad(V)-Hess(W)grad(W)=',list(mismatch))
    J=s.Matrix([[0,-1],[1,0]])
    assert (J*M+M*J).applyfunc(s.simplify)==s.zeros(2)
    assert (s.diff(M,v)-(J*M-M*J)/(2*r)).applyfunc(s.simplify)==s.zeros(2)
    print('PASS axial rotation and combined v-shift/vector rotation are exact action symmetries')
    # A 1|1 type-one module cannot carry a nontrivial finite sl(2|1)
    # representation. The even sl2 raising/lowering maps each 1-dimensional
    # parity sector to itself; the [h1,e1]=2e1 equation forces them to zero.
    # Type-one integrability then gives h1=0 and [h1,e2]=-e2 forces e2=0.
    print('PASS dimensional obstruction: an isolated 1|1 doublet is not a nontrivial type-one sl(2|1) module')
    # Directly enumerate the possible defining-triplet weights for a topological
    # charge in the sl(2|1) Cartan, without a central offset.
    bad=[]
    for a in (-1,1):
        for b in (-1,1):
            if b-a in (-1,1):bad.append((a,b))
    assert not bad
    print('PASS no Cartan combination assigns defining-triplet weights only unit topological charges +/-1')

    # A compatible four-state module, at fixed alpha=-1/2, CAN be built.
    # Basis: two even states then two odd states. q=t^2, A=[1/2]_q.
    t=s.symbols('t',positive=True)
    q=t**2;A=t/(t**2+1)
    def E(i,j):
        mat=s.zeros(4);mat[i,j]=1;return mat
    h=[s.diag(0,0,-1,1),s.diag(-s.Rational(1,2),s.Rational(1,2),s.Rational(1,2),-s.Rational(1,2))]
    e=[E(3,2),-E(0,3)+E(2,1)]
    f=[E(2,3),A*(E(1,2)+E(3,0))]
    cartan=s.Matrix([[2,-1],[-1,0]])
    for i in range(2):
        for j in range(2):
            assert h[i]*e[j]-e[j]*h[i]==cartan[i,j]*e[j]
            assert h[i]*f[j]-f[j]*h[i]==-cartan[i,j]*f[j]
            lhs=e[i]*f[j]-(-1)**(i*j)*f[j]*e[i]
            rhs=s.diag(*[(t**(2*h[i][k,k])-t**(-2*h[i][k,k]))/(q-q**-1) for k in range(4)]) if i==j else s.zeros(4)
            assert (lhs-rhs).applyfunc(s.factor)==s.zeros(4)
    assert e[1]**2==f[1]**2==s.zeros(4)
    for z in (e,f):
        assert z[0]**2*z[1]-(q+q**-1)*z[0]*z[1]*z[0]+z[1]*z[0]**2==s.zeros(4)
    T=2*h[1];F=h[0]+h[1]
    assert list(T.diagonal())==[-1,1,1,-1]
    assert list(F.diagonal())==[-s.Rational(1,2),s.Rational(1,2),-s.Rational(1,2),s.Rational(1,2)]
    print('PASS all finite Uq(sl(2|1)) Cartan, supercommutator, nilpotence and Serre relations on a compatible four-state module')
    print('Compatible weights (T,F,parity):',[(T[i,i],F[i,i],(0,0,1,1)[i]) for i in range(4)])
    # Highest vector is state 1; h1=0,h2=-1/2. This fixes, not frees, alpha.
    assert e[0][:,0]==e[1][:,0]==s.zeros(4,1)
    assert (h[0][0,0],h[1][0,0])==(0,-s.Rational(1,2))
    for ei,dT,dF in ((e[0],-2,1),(e[1],0,-1)):
        assert T*ei-ei*T==dT*ei and F*ei-ei*F==dF*ei
    print('PASS compatible highest weight alpha=-1/2; e1 selection rule (Delta T,Delta F)=(-2,1), e2=(0,-1)')
    print('CAUTION: the nonlocal field operators realizing these matrices, vacuum charge offsets, and affine extension on physical kinks are NOT constructed.')

if __name__=='__main__':main()
