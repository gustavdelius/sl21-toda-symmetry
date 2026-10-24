"""Fresh Wick/Taylor calculation. No inherited scripts imported."""
import sympy as S
import mpmath as mp

t,z,a,r,v,q,eps=S.symbols('t z a r v q eps', nonzero=True)
rho=S.symbols('rho')
b=S.symbols('b', nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I])
a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I])
T=a0+a2; d=-T/b**2
count=0
def check(name, cond):
    global count
    count+=1
    print(('PASS' if bool(cond) else 'FAIL')+' '+name)
    if not cond: raise AssertionError(name)
check('bilinear root and disorder products', [S.simplify(x) for x in [a0.dot(a0)+1,a2.dot(a2)+1,a0.dot(a2)-(1-b*b),a0.dot(d)-1,a2.dot(d)-1,d.dot(d)+2/b**2]]==[0]*6)
check('K2/Y0 chiral exponent t-2',S.simplify(-(d+a2).dot(a0)-(b*b-2))==0)
check('U/Y0 full exponent 2t-4',S.simplify(-2*(d+a2).dot(a0)-(2*b*b-4))==0)
check('source pair exponent 2t-2',S.simplify(-2*a2.dot(a0)-(2*b*b-2))==0)
f=(1/(z-r)+(t-1)/z)/((z-r)*z)
h=(a-r)*a*(1/(z-a)+1/(z-r)+(t-1)/z)/((z-r)*z)
# h suppresses only the spectator factor (a-z)^(-rho).
for n in range(6):
    fn=S.diff(f,r,n).subs(r,0)/S.factorial(n)
    hn=S.diff(h,r,n).subs(r,0)/S.factorial(n)
    check(f'bare chiral coefficient n={n}',S.simplify(fn-(n+t)/z**(n+3))==0)
    target=a*a/z**3*(t-z/(a-z)) if n==0 else a*(a-z)*(n+t)/z**(n+3)
    check(f'actual A-dressed chiral coefficient n={n}',S.simplify(hn-target)==0)
closed=q*((1+t)-t*q)/(1-q)**2
check('exact source-pair remainder sum',S.simplify(S.diff(closed,q)-(1+t+(1-t)*q)/(1-q)**3)==0)
check('leading remainder coefficient nonzero generic',S.simplify(S.limit(closed/q,q,0)-(1+t))==0)
# Minus action insertion yields the displayed negative cutoff representatives.
lam0,lam2,al=S.symbols('lambda0 lambda2 alpha')
lead=-lam0/(t-1)*eps**(2*t-2)
desc=-lam0/t*eps**(2*t)
check('canonical endpoint leading subtraction',S.simplify(t*t*lam2*lead+t*t*lam0*lam2/(t-1)*eps**(2*t-2))==0)
check('canonical endpoint descendant subtraction',S.simplify(t*t*lam2*desc+t*lam0*lam2*eps**(2*t))==0)
check('finite matching leading subtraction',S.simplify(t*(t-1)*lam2*lead+t*lam0*lam2*eps**(2*t-2))==0)
check('finite matching descendant subtraction',S.simplify(t*(t-1)*lam2*desc+(t-1)*lam0*lam2*eps**(2*t))==0)
check('tail improvement coefficient',S.simplify((-t*lam2)*(-lam0/t)*eps**(2*t)-lam0*lam2*eps**(2*t))==0)
labels={'A':(-1,0,0,0),'B':(1,0,0,0),'b':(0,1,0,0),'0':(0,0,1,0),'2':(0,0,0,1),'P':(0,0,1,1),'K':(1,0,0,1),'k':(0,1,0,1),'L':(1,0,1,1),'l':(0,1,1,1)}
def cocycle_word(word):
    sign=0; phase=0
    for i,nm in enumerate(word):
        n=labels[nm]
        for mm in word[i+1:]:
            m=labels[mm]
            sign+=(n[0]+n[1])*(m[2]+m[3])
            phase+=n[0]*m[1]-n[1]*m[0]
    return (-1 if sign%2 else 1),phase
for word,expected in [('ABb02',(1,0)),('ABbP',(1,0)),('ABk0',(-1,0)),('AKb0',(1,0)),('ABl',(1,0)),('ALb',(1,0))]:
    check('reference-order cocycle '+word,cocycle_word(word)==expected)

mp.mp.dps=40
tv=-mp.sqrt(2)/3
av=mp.mpf(3)
# Positive heights compatible with one X14 word; bar current remains at height zero.
# w2=0, currents approach x=M epsilon; A is held away from this proper face.
M=mp.mpf(12); c=mp.sqrt(5)
print('t =',mp.nstr(tv,20),'c =',mp.nstr(c,20),'M =',M)
def exact(e,zv,vv):
    Q=e*e/(zv*vv)
    return mp.pi*e**(2*tv)/(zv*vv)**3*Q*((1+tv)-tv*Q)/(1-Q)**2
def series_integral(e,zv,vv):
    # Integrate each angular mode independently, retaining fifty nonzero modes.
    return mp.fsum([mp.pi*(n+tv)*e**(2*tv+2*n)/(zv*vv)**(3+n) for n in range(1,51)])
for e in [mp.mpf('0.08'),mp.mpf('0.04'),mp.mpf('0.02'),mp.mpf('0.01')]:
    zv=e*(M+2j*c); vv=e*M
    E=exact(e,zv,vv)
    angular=mp.fsum([((1/(zv-e*mp.e**(2j*mp.pi*k/2048))+(tv-1)/zv)/((zv-e*mp.e**(2j*mp.pi*k/2048))*zv))*((1/(vv-e*mp.e**(-2j*mp.pi*k/2048))+(tv-1)/vv)/((vv-e*mp.e**(-2j*mp.pi*k/2048))*vv)) for k in range(2048)])/2048
    Q=e*e/(zv*vv)
    expected=mp.fsum([(n+tv)**2*Q**n for n in range(50)])/(zv*vv)**3
    check('angular quadrature epsilon='+str(e),abs((angular-expected)/expected)<mp.mpf('1e-30'))
    check('integrated mode sum epsilon='+str(e),abs((E-series_integral(e,zv,vv))/E)<mp.mpf('1e-30'))
    print('DIAGNOSTIC epsilon',e,'absolute omitted disk',mp.nstr(abs(E),15),'scaled coefficient',mp.nstr(E*e**(6-2*tv),15))
e1=mp.mpf('.02');e2=e1/2
R1=exact(e1, e1*(M+2j*c),e1*M)
R2=exact(e2, e2*(M+2j*c),e2*M)
check('nonuniform intersecting-face scaling',abs(R2/R1-2**(6-2*tv))<mp.mpf('1e-30'))
fixed1=exact(e1,2+0.5j,2)
fixed2=exact(e2,2+0.5j,2)
print('DIAGNOSTIC fixed-separation decay ratio',mp.nstr(abs(fixed2/fixed1),15),'limit',mp.nstr(2**(-2*tv-2),15))
check('fixed-separation remainder decays',abs(fixed2)<abs(fixed1))
print('Completed',count,'substantive checks; this is not an evaluated four-word defect.')
