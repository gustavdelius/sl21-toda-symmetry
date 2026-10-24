"""Independent defining/dual Uq(sl(2|1)^(1)) intertwiners.

Column-vector convention, parity (0,0,1), q=t**2, inverse evaluation z.
The matrix x below is the ratio z_left/z_right; e0(z)=z^-1 E31.
All tensor products are realised as ordinary matrices with Koszul signs explicit.
This is referee scratch, not book code. No imports from the author's checks.
"""
from pathlib import Path
import pickle
import sympy as s

t, x = s.symbols('t x', nonzero=True)
q = t**2
p = (0, 0, 1)
J = s.diag(1, 1, -1)
I = s.eye(3)

def E(i, j):
    a = s.zeros(3); a[i,j] = 1
    return a

def st(a):
    # (E_ij)^st = (-1)^(p_i*(p_i+p_j)) E_ji.
    return s.Matrix(3,3,lambda i,j: (-1)**(p[j]*(p[j]+p[i]))*a[j,i])

def gt(a,b, pb=0):
    return s.kronecker_product(a*(J if pb else I), b)

P = s.zeros(9)
for i in range(3):
    for j in range(3):
        P[3*j+i,3*i+j] = (-1)**(p[i]*p[j])

def rep(z, dual=False):
    h = [s.diag(-1,0,-1),s.diag(1,-1,0),s.diag(0,1,1)]
    e = [E(2,0)/z,E(0,1),E(1,2)]
    f = [-z*E(0,2),E(1,0),E(2,1)]
    if dual:
        # rho*(a)=rho(S(a))^st. S(e)=-K^-1 e K, likewise f.
        es=[]; fs=[]
        for hi,ei,fi in zip(h,e,f):
            K=s.diag(*[t**hi[i,i] for i in range(3)])
            es.append(st(-K.inv()*ei*K))
            fs.append(st(-K.inv()*fi*K))
        e,f,h=es,fs,[-a for a in h]
    return {'h':h,'e':e,'f':f}

def coproduct(a,b,key,i):
    if key=='h':
        return gt(a[key][i],I)+gt(I,b[key][i])
    pa=(1,0,1)[i]
    ka=s.diag(*[t**a['h'][i][k,k] for k in range(3)])
    kb=s.diag(*[t**b['h'][i][k,k] for k in range(3)])
    return gt(a[key][i],kb.inv())+gt(ka,b[key][i],pa)

def intertwiner(da,db):
    a,b=rep(x,da),rep(1,db)
    weights=[tuple(a['h'][k][i,i]+b['h'][k][j,j] for k in (1,2))
             for i in range(3) for j in range(3)]
    positions=[(i,j) for i in range(9) for j in range(9) if weights[i]==weights[j]]
    vs=s.symbols('r:'+str(len(positions)))
    R=s.zeros(9)
    for ij,v in zip(positions,vs): R[ij[0],ij[1]]=v
    eq=[]
    for key in ('e','f'):
        for k in range(3):
            D=coproduct(a,b,key,k)
            Dop=P*coproduct(b,a,key,k)*P
            eq.extend(v for v in R*D-Dop*R if v!=0)
    norm=0 if da==db else 1
    eq.append(R[norm,norm]-1)
    mat,rhs=s.linear_eq_to_matrix(eq,vs)
    sol=s.linsolve((mat,rhs),vs)
    tup=next(iter(sol))
    out=R.subs(dict(zip(vs,tup))).applyfunc(s.factor)
    assert not (set(vs)&out.free_symbols)
    return out

def paper_check():
    # Independently solved above; transcription used only for comparison.
    r=s.zeros(9); r[0,0]=r[4,4]=1
    r[8,8]=(1-x*q**2)/(x-q**2)
    A=(1-q**-2)/(1-x*q**-2)
    B=q**-1*(1-x)/(1-x*q**-2)
    for j in (1,2,5): r[j,j]=A
    for j in (3,6,7): r[j,j]=x*A
    for i,j in ((1,3),(2,6),(5,7),(3,1),(6,2),(7,5)): r[i,j]=B
    return r

def main():
    folder=Path(__file__).parent
    mats={}
    for da,db in ((False,False),(True,False),(False,True),(True,True)):
        name=('b' if da else 'v')+('b' if db else 'v')
        print('Solving',name,flush=True)
        R=intertwiner(da,db); mats[name]=R
        print(name,[(i,j,R[i,j]) for i in range(9) for j in range(9) if R[i,j]!=0],flush=True)
    delta=(P*mats['vv']-paper_check()).applyfunc(s.factor)
    assert delta==s.zeros(9),delta
    print('PASS independent affine intertwiner equals (4.11) after graded permutation')
    with (folder/'matrices.pkl').open('wb') as f: pickle.dump(mats,f)

if __name__=='__main__': main()
