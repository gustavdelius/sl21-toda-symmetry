"""Independent dual, mixed crossing and matrix-unitarity calculations."""
from algebra import *

def clean(a): return a.applyfunc(s.factor)
def st1(a):
    # Translate the *graded tensor* supertranspose into ordinary matrix storage.
    # A term E_ik tensor_g E_jl stores (-1)^((p_j+p_l)*p_k).
    out=s.zeros(9)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    out[3*k+j,3*i+l]=(-1)**((p[j]+p[l])*(p[k]+p[i])+p[i]*(p[i]+p[k]))*a[3*i+j,3*k+l]
    return out

def dual_of(a):
    b={'h':[-h for h in a['h']],'e':[],'f':[]}
    for k in range(3):
        K=s.diag(*[t**a['h'][k][i,i] for i in range(3)])
        for key in ('e','f'): b[key].append(clean(st(-K.inv()*a[key][k]*K)))
    return b

def load():
    with (Path(__file__).parent/'matrices.pkl').open('rb') as f: return pickle.load(f)

def main():
    mats=load()
    double=dual_of(rep(x,True)); D=s.diag(q**2,1,-1)
    for key in ('h','e','f'):
        for k in range(3):
            assert clean(D*double[key][k]*D.inv()-rep(q**2*x)[key][k])==s.zeros(3)
    print('PASS double dual is D^-1 V(q^2*z) D, D=diag(q^2,1,-1)')
    for name,R in mats.items():
        prod=clean(R*P*mats[name[::-1]].subs(x,1/x)*P)
        scalar=s.factor(prod[0,0]); assert clean(prod-scalar*s.eye(9))==s.zeros(9)
        print('PASS matrix braiding unitarity',name,'scalar=',scalar)
    # (S tensor id) R = R^-1, hence R_bv proportional to st1(R_vv^-1).
    rhs=clean(st1(mats['vv'].inv()))
    fac=s.factor(mats['bv'][1,1]/rhs[1,1])
    assert clean(mats['bv']-fac*rhs)==s.zeros(9)
    print('PASS R_bv(x) =',fac,'* st1(R_vv(x)^-1)')
    # A second application also checks the convention for the nontrivial st1.
    rhs=clean(st1(mats['bv'].inv()))
    DD=s.kronecker_product(D,s.eye(3))
    lhs=clean(DD*rhs*DD.inv())
    shifted=mats['vv'].subs(x,q**2*x)
    fac2=s.factor(lhs[0,0]/shifted[0,0])
    assert clean(lhs-fac2*shifted)==s.zeros(9)
    print('PASS second dual crossing scalar=',fac2)
    rhs=clean(st1(mats['vb'].inv()))
    fac3=s.factor(rhs[0,0]/mats['bb'][0,0])
    assert clean(rhs-fac3*mats['bb'])==s.zeros(9)
    print('PASS st1(R_vb^-1) =',fac3,'* R_bb')
    rhs=clean(DD*st1(mats['bb'].inv())*DD.inv())
    fac4=s.factor(rhs[1,1]/mats['vb'].subs(x,q**2*x)[1,1])
    assert clean(rhs-fac4*mats['vb'].subs(x,q**2*x))==s.zeros(9)
    print('PASS st1(R_bb^-1), double-dual identified =',fac4,'* R_vb(q^2*x)')
    # Trivial and adjoint components of V* tensor V at finite-algebra level.
    a,b=rep(1,True),rep(1)
    sing=s.Matrix([1,0,0,0,1,0,0,0,1])
    ds=[coproduct(a,b,key,k) for key in ('h','e','f') for k in (1,2)]
    basis=s.Matrix.vstack(*ds).nullspace()
    print('finite singlet basis=',[list(v) for v in basis])
    assert len(basis)==1
    # Characteristic polynomial; all exact rational functions.
    for name in ('bv','vb'):
        print(name,'unpermuted eigenvalues=',mats[name].eigenvals())

if __name__=='__main__': main()
