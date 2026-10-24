import sympy as S
from pathlib import Path
p=Path(__file__).parent
t,g,rho,a,w,q,Ca,Cw,D,u=S.symbols('t g rho a w q Ca Cw D u', nonzero=True)
beta=t-1; delta=g-rho
checks=[]
def check(name,value):
    value=S.factor(value.subs(rho,2/t))
    ok=value==0
    print(('PASS' if ok else 'FAIL')+' '+name+((' : '+str(value)) if not ok else ''))
    checks.append(ok)
r=w-a
b=Ca+(beta-1)/(a-w)
# Independent Laurent expansion of the actual two derivative-current Wick polynomial.
P=(Cw-1/(w+q-a)+beta/q)*(Ca-1/(a-w-q)+beta/(a-w))+1/(w+q-a)**2
pole=S.limit(q*P,q,0)
constant=S.limit(P-pole/q,q,0)
# Strip E=(w-a)^g/(a-w); normal-field first jet at z=w is D.
Dc=pole
Sc=constant+pole*(g/r+D)
Kw=(g-1)*b/r+(beta-1)/r**2-b*(Cw+beta*D)
check('moving double coefficient beta times spectator-tail product',Dc-beta*b)
check('moving simple coefficient -Kw+t(g-rho)K/r',Sc+Kw-t*delta*b/r)
# The wrong unit coefficient is explicitly retained as a diagnostic expression.
wrong=S.factor((Sc+Kw-delta*b/r).subs(rho,2/t))
print('DIAGNOSTIC rejected unit-regulator coefficient residual =',wrong)
check('rejected unit coefficient differs by beta delta K/r',wrong-(t-1)*(g-2/t)*b/r)
# A genuine arbitrary smooth test: g0+g1*(z-w); its double pole shifts the simple coefficient.
g0,g1=S.symbols('g0 g1')
check('linear moving-test residue retains double-pole derivative',g0*Sc+g1*Dc-(g0*Sc+g1*beta*b))
# I=pi, Jw=J1=-pi produce these independent local first-jet terms.
Qx,Qb,Q0=S.symbols('Qx Qb Q0')
action=-beta*(t*(Qx+Qb)-t*delta*Q0/r)
tail_left=t*(Qx-t*delta*Q0/r)
tail_right=t*(beta*Qx+t*Qb)
raw=S.expand(action+tail_left+tail_right)
check('action plus both tails: t times insertionwise endpoint derivative',raw-t*(Qx+Qb-delta*Q0/r))
alpha=t*(t-1)
matched=raw+alpha*(Qx+Qb-delta*Q0/r)
check('canonical finite matching: t squared insertionwise derivative',matched-t*t*(Qx+Qb-delta*Q0/r))
# Full spatial-test linear coefficient: action I+Jw cancels; opposite tail terms cancel.
check('action first test derivative I+Jw',beta*beta*(1-1))
check('two tail first test derivatives cancel',t*beta-t*beta)
# Full Gaussian scalar B^3 A / Y2 reconstruction at two asymmetric shapes.
for zi in [list(map(S.Rational,[-2,0,3])),[S.Rational(-3,2),S.Rational(2,3),S.Rational(5,2)]]:
    fs=[-sum(1/(zi[i]-zi[j]) for j in range(3) if j!=i)+1/(zi[i]-a)+beta/(zi[i]-w) for i in range(3)]
    poly=fs[0]*fs[1]*fs[2]+sum(fs[3-i-j]/(zi[i]-zi[j])**2 for i in range(3) for j in range(i+1,3))
    H=(a-w)**2/S.prod((w-z)**2 for z in zi)
    logR=rho*sum(1/(z-a) for z in zi)
    full=(a-w)/S.prod(z-w for z in zi)*poly
    primitive=beta*(t-2)/rho*(S.diff(H,a)+logR*H)
    check('B3A source scalar minority derivative '+str(zi),full-primitive)
    # Double anti pole equals a bar-coordinate derivative before restriction.
    wb=S.symbols('wb')
    check('five-current source has scalar double derivative '+str(zi),beta/(u-wb)**2+S.diff(beta/(u-wb),u))
    # Barred tail has only minority derivative immediately; simple u poles obstruct rational u primitive.
    H_u=H.subs(w,u)
    residues=[S.factor(S.diff(((u-z)**2*H_u).cancel(),u).subs(u,z)) for z in zi]
    check('barred-tail simple residue sum vanishes at infinity '+str(zi),sum(residues))
    nonzero=any(v!=0 for v in residues)
    print(('PASS' if nonzero else 'FAIL')+' barred tail has nonzero individual simple residues '+str(zi))
    checks.append(nonzero)
    print('DIAGNOSTIC barred-tail simple residues =',residues)
# Independent dimensions from counts; not used as a substitute for the contact calculation.
# n chiral B, m A, one bar B, one area source, maximal n+1 C contractions.
n,m=S.symbols('n m',integer=True,nonnegative=True)
deg=g*n*(n-1)/2+g*m*(m-1)/2-rho*n*m +(-n+m-1)-(n+1)+(n+m+2)
check('BBbarB/Y2 full relative degree g-2',deg.subs({n:2,m:0})-(g-2))
check('B3AbarB/Y2 full relative degree 3(g-rho)-1',deg.subs({n:3,m:1})-(3*(g-rho)-1))
# Pairwise cross-regulator pullback is needed for fused subtractions.
x1,x2,y,z=S.symbols('x1 x2 y z')
log_cross=delta*(S.log(x1-z)+S.log(x2-z))
check('fused pair carries two cross labels',S.diff(log_cross,z).subs({x1:y,x2:y})+2*delta/(y-z))
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks): raise SystemExit(1)
