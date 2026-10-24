"""A crossing-compatible algebraic six-state scattering candidate and first fusion.

This does not assert Hermitian analyticity, positive-metric unitarity, or a Toda
identification. q=exp(i*pi*kappa), X=exp(kappa*theta), 0<kappa<1 is a
pole-free branch. Barred particles use the left dual at evaluation z/q.
"""
from crossing import *
import numpy as np

def scalar(name,z):
    if name in ('vv','bb'): return s.S.One
    if name=='bv': return q*(z-1)/(q**2*z-1)
    if name=='vb': return (z-q**2)/(q*(z-1))
    raise KeyError(name)

def raw_candidate():
    return {n:clean(scalar(n,x)*r) for n,r in load().items()}

def physical_candidate():
    mats=raw_candidate()
    return {n:clean(r.subs(x,x*(q if n=='vb' else 1/q if n=='bv' else 1)))
            for n,r in mats.items()}

def main():
    raw=raw_candidate(); phys=physical_candidate(); DD=s.kronecker_product(s.diag(q**2,1,-1),I)
    for n,r in raw.items():
        assert clean(r*P*raw[n[::-1]].subs(x,1/x)*P)==s.eye(9)
        target=raw[('b' if n[0]=='v' else 'v')+n[1]]
        crossed=clean(st1(r.inv()))
        if n[0]=='b':
            crossed=clean(DD*crossed*DD.inv());target=target.subs(x,q**2*x)
        assert clean(crossed-target)==s.zeros(9)
        print('PASS exact raw braiding unitarity and left-dual crossing',n,flush=True)
    for n,r in phys.items():
        assert clean(r*P*phys[n[::-1]].subs(x,1/x)*P)==s.eye(9)
        print('PASS exact physical braiding unitarity',n,flush=True)
        den=s.factor(s.lcm([s.denom(s.cancel(z)) for z in r]))
        print('physical common denominator',n,den,flush=True)
    # All eight species assignments in the graded YBE. All reps have same parity.
    P23=np.kron(np.eye(3),np.array(P,float))
    f={n:s.lambdify((t,x),r,'numpy') for n,r in phys.items()}
    worst=0.
    for kap in (np.sqrt(2)/3,np.sqrt(3)/2):
        tn=np.exp(0.5j*np.pi*kap)
        for a in 'vb':
            for b in 'vb':
                for c in 'vb':
                    r12=np.kron(f[a+b](tn,np.exp(kap*.47)),np.eye(3))
                    r23=np.kron(np.eye(3),f[b+c](tn,np.exp(kap*-.31)))
                    r13=P23@np.kron(f[a+c](tn,np.exp(kap*.16)),np.eye(3))@P23
                    err=np.max(np.abs(r12@r13@r23-r23@r13@r12));worst=max(worst,err)
    assert worst<1e-11,worst
    print('PASS all 8 mixed graded YBE assignments, 2 irrational kappas; max residual',worst,flush=True)
    # Grade/rapidity changes cannot be inferred from the symmetry alone.
    tn=np.exp(0.5j*np.pi*np.sqrt(2)/3);xn=np.exp(np.sqrt(2)/3*.47)
    B=np.array(P,float)@f['vv'](tn,xn)
    sv=np.linalg.svd(B,compute_uv=False)
    print('Homogeneous vv singular values at real rapidity:',sv,flush=True)
    assert np.max(np.abs(sv-1))>1e-3
    print('PASS diagnostic: braiding unitarity is not positive-metric unitarity',flush=True)
    # Exact rank-four residue at the raw vv pole x=q^2.
    R=load()['vv']; check=P*R
    residue=clean(((x-q**2)*check).applyfunc(s.cancel).subs(x,q**2))
    assert residue.rank()==4
    assert clean(check*check.subs(x,1/x))==s.eye(9)
    print('PASS vv pole residue has rank 4 (not antifundamental rank 3)',flush=True)
    # First fusion in inverse evaluation convention. Image is in exchanged order.
    B=s.zeros(9,4)
    for k,(i,j) in enumerate(((0,1),(0,2),(1,2))):
        B[3*i+j,k]=1; B[3*j+i,k]=-q
    B[8,3]=1
    L=s.zeros(4,9)
    for k,j in enumerate((1,2,5,8)): L[k,j]=1
    assert L*B==s.eye(4)
    # Highest vector in the rank-four channel has (h1,h2)=(0,1), so it
    # is the alpha=1 member of (0,0|alpha) on restriction to sl(2|1).
    highest=B[:,0]
    for j in (1,2):
        assert clean(coproduct(rep(1),rep(1),'e',j)*highest)==s.zeros(9,1)
    assert clean(coproduct(rep(1),rep(1),'h',1)*highest)==s.zeros(9,1)
    assert clean(coproduct(rep(1),rep(1),'h',2)*highest-highest)==s.zeros(9,1)
    print('PASS rank-four fusion highest weight (h1,h2)=(0,1): alpha=1, no free continuous label',flush=True)
    B3=s.kronecker_product(B,I); L3=s.kronecker_product(L,I)
    P23s=s.kronecker_product(I,P)
    R13=P23s*s.kronecker_product(R.subs(x,-x/q),I)*P23s
    R23=s.kronecker_product(I,R.subs(x,-q*x))
    full=R13*R23
    fused=clean(L3*full*B3)
    leak=clean(full*B3-B3*fused)
    assert leak==s.zeros(27,12)
    print('PASS first fused 4x3 amplitude preserves the rank-four bound-state space',flush=True)
    den=s.factor(s.lcm([s.denom(s.cancel(z)) for z in fused]))
    print('fused common denominator:',den,flush=True)
    for pole in (-q,-q**3):
        res=clean(((x-pole)*fused).applyfunc(s.cancel).subs(x,pole))
        print('fused residue at',pole,'rank',res.rank(),flush=True)
    # Scatter the same fused state against an antifundamental, including its
    # crossing scalar and the physical dual evaluation offset.
    C=phys['vb']
    C13=P23s*s.kronecker_product(C.subs(x,-x/q),I)*P23s
    C23=s.kronecker_product(I,C.subs(x,-q*x))
    fullbar=C13*C23
    fusedbar=clean(L3*fullbar*B3)
    assert clean(fullbar*B3-B3*fusedbar)==s.zeros(27,12)
    print('PASS fused 4 x antifundamental preserves the bound-state space',flush=True)
    print('fused 4 x antifundamental denominator:',s.factor(s.lcm([s.denom(s.cancel(z)) for z in fusedbar])),flush=True)
    for pole in (-s.S.One,-1/q**2):
        res=clean(((x-pole)*fusedbar).applyfunc(s.cancel).subs(x,pole))
        print('fused 4 x antifundamental residue at',pole,'rank',res.rank(),flush=True)
    with (Path(__file__).parent/'candidate.pkl').open('wb') as file:
        pickle.dump({'raw':raw,'physical':phys,'fused43':fused,'fused4b3':fusedbar},file)

if __name__=='__main__': main()
