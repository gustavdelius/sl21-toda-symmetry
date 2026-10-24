"""An exact Gaussian probe with the two physical bosons held fixed.

This is NOT the expectation value of the full renormalized endpoint U.
It tests only the invalid inference from classical psi=0 to zero bilinear VEV.
"""
import sympy as S
import mpmath as mp

n=0
def check(name,v):
    global n
    if not bool(v):
        print('FAIL',name)
        raise AssertionError(name)
    n+=1
    print('PASS',name)
def eq(e): return S.simplify(e)==0

p,pb,m0,m2,p2,Lambda,M,mu=S.symbols('p pb m0 m2 p2 Lambda M mu',nonzero=True)
A=S.Matrix([[0,pb,S.I*m0,0],[pb,0,0,S.I*m2],[-S.I*m0,0,0,p],[0,-S.I*m2,p,0]])
check('fermionic kernel has A(p)=-A(-p)^T',A+A.subs({p:-p,pb:-pb}).T==S.zeros(4))
check('mass pole from determinant',eq(A.det()-(p*pb-m0*m2)**2))
cov=S.pi*A.inv()
check('covariance normalizes the free chiral OPE',eq(cov[0,1].subs({m0:0,m2:0})-S.pi/pb))
entry=S.factor(cov[1,3])
euclidean=entry.subs(pb,-p2/(4*p))
check('chi2-bar chi2 massive propagator',eq(euclidean-4*S.pi*S.I*m0/(p2+4*m0*m2)))
check('i chi2 bar chi2 has negative momentum kernel',eq(S.I*euclidean+4*S.pi*m0/(p2+4*m0*m2)))
cutoff=-m0*S.log((Lambda**2+M**2)/M**2)
check('momentum radial primitive',eq(S.diff(cutoff,Lambda)+2*m0*Lambda/(Lambda**2+M**2)))
check('subtracted finite part',eq(S.limit(cutoff+m0*S.log(Lambda**2/mu**2),Lambda,S.oo)-m0*S.log(M**2/mu**2)))

mp.mp.dps=40
mass=mp.sqrt(2); coeff=mp.e/3
for r in [mp.mpf('.3'),mp.mpf('.7'),mp.sqrt(3)]:
    # Schwinger-parameter Fourier transform, independently evaluated.
    green=mp.quad(lambda x: mp.exp(-mass**2*x-r*r/(4*x))/(4*mp.pi*x),[0,1,mp.inf])
    expected=mp.besselk(0,mass*r)/(2*mp.pi)
    check(f'Fourier transform at r={mp.nstr(r,7)}',abs(green-expected)<mp.mpf('1e-32'))
    check(f'bilinear Bessel coefficient at r={mp.nstr(r,7)}',abs(-4*mp.pi*coeff*green+2*coeff*mp.besselk(0,mass*r))<mp.mpf('1e-31'))
reference=mp.sqrt(3)
target=2*coeff*(mp.log(mass/reference)+mp.euler-mp.log(2))
errs=[]
for r in [mp.mpf('1e-3'),mp.mpf('1e-4'),mp.mpf('1e-5')]:
    finite=-2*coeff*mp.besselk(0,mass*r)-2*coeff*mp.log(reference*r)
    errs.append(abs(finite-target))
check('coordinate finite part contains the stated scheme constant',errs[-1]<mp.mpf('2e-9') and errs[-1]<errs[0]/1000)
print('Coordinate finite-part errors:',', '.join(mp.nstr(x,10) for x in errs))
print(f'{n} PASS, 0 FAIL')
