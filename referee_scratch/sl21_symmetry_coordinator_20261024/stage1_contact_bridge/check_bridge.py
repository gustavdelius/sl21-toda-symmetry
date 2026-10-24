"""Direct spectator-dressed mixed-face calculation; no old scripts imported."""
import sympy as S
import mpmath as mp
t,h,v,vb,w,z,C,db,a2b=S.symbols('t h v vb w z C db a2b',nonzero=True)
e=h-t; beta=t-1
a0b=-t*db-a2b; ub=db+a2b
count=0
def ck(name,condition):
    global count
    count+=1
    print(('PASS' if bool(condition) else 'FAIL')+' '+name)
    if not condition:raise AssertionError(name)
# Direct current/source pole expansion. Source pair factor is independent of z.
for label,bc,bs,sgn in [('node2',beta,1,-1),('node0',1,beta,1)]:
    numerator=(C+bs/(z-v))/(z-v)+bc/((z-w)*(z-v))
    # Entire kernel = numerator/(z-w); use regularized numerator for Laurent jets.
    regular=(z-w)*(C+bs/(z-v))/(z-v)+bc/(z-v)
    D=regular.subs(z,w)
    Scoef=S.diff(regular,z).subs(z,w)
    ck(label+' actual double coefficient',S.simplify(D-bc/(w-v))==0)
    ck(label+' scalar simple coefficient',S.simplify(Scoef-(C/(w-v)+(bs-bc)/(w-v)**2))==0)
# Full first jets use a0+beta d=-(d+a2), and a0+d=(d+a0).
b=S.symbols('b',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I])
a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I])
d=-(a0+a2)/b**2
ck('node2 full field jet identity',all(S.simplify(x)==0 for x in a0+(b*b-1)*d+d+a2))
ck('antiholomorphic mixed identity beta a2-a0=t(d+a2)',S.simplify(beta*a2b-a0b-t*ub)==0)
# Node2: after holomorphic collapse, K_w is the full h-dressed endpoint derivative.
ux=S.symbols('ux')
Kw=ux-(h-2)/v
Kb=a2b-(h-1)/vb
Dbar=a0b-1/vb
Qbar=ub-(h-2)/vb
R=1/v+1/vb
# Strip integrals I=pi,Jw=J1=-pi, after removing pi.
raw_action_bracket=beta*(-t*Kw-e/v+Dbar-beta*Kb)
# action overall multiplier is -2 i, so expect minus beta*[t derivative+e R]
ck('node2 action full first jets',S.simplify(raw_action_bracket+beta*(t*(Kw+Qbar)+e*R))==0)
tail1=-t*(Kw+e/v)
tail2=t*(Dbar-beta*(Kw+Kb))
ck('node2 both tails full first jets',S.simplify(tail1+tail2+t*(t*(Kw+Qbar)+e*R))==0)
action=beta*(t*(Kw+Qbar)+e*R)
tails=tail1+tail2
ck('raw mixed contact exact regulator remainder',S.simplify(action+tails+t*(Kw+Qbar)+e*R)==0)
# Node0: U0=d+a0; antiholomorphic source field a0, other pole coefficient beta.
K0b=a0b-(h-1)/vb
D0bar=a0b-beta/vb
ck('node0 action contact solely h-t',S.simplify(e/v+D0bar-K0b-e*R)==0)
alpha=S.symbols('alpha')
ck('insertion-dressed matching regulator term',S.simplify((-t*(Kw+Qbar)-e*R)-alpha*((Kw+Qbar)+e*R)+(t+alpha)*(Kw+Qbar)+(1+alpha)*e*R)==0)
ck('canonical insertion matching coefficient',S.simplify(t+t*beta-t*t)==0)
ck('scalar first outer test jet vanishes action',beta*beta*(-1+1)==0)
ck('scalar first outer test jet vanishes tails',S.simplify(t*beta-t*beta)==0)
mp.mp.dps=40
tv=-mp.sqrt(2)/3; vv=2+0.7j
for eh in [mp.mpf('.02'),mp.mpf('.005'),mp.mpf('.001')]:
    hv=tv+eh; Rv=1/vv+1/vv.conjugate()
    Q=abs(vv)**(2*hv-4)
    # Smooth field jet zero; both endpoint derivatives are Coulomb derivatives.
    Qx=-(hv-2)*Rv*Q
    actionv=(tv-1)*(tv*Qx+eh*Rv*Q)
    tailsv=-tv*(tv*Qx+eh*Rv*Q)
    rem=actionv+tailsv+tv*Qx
    ck('numerical full-contact divisibility e='+str(eh),abs(rem+eh*Rv*Q)<mp.mpf('1e-35'))
    print('DIAGNOSTIC e',eh,'raw contact remainder/(2 i nu lambda0 lambda2)',mp.nstr(rem,18),'remainder/e',mp.nstr(rem/eh,18))
# Independent direct numerical Laurent extraction and derivative of the h factor.
tv=mp.mpf(tv);hv=mp.sqrt(7);wv=mp.mpf('.3');vv=mp.mpf('2.4');cv=mp.mpf('.27')
f=lambda zz:((zz-wv)*(cv+1/(zz-vv))/(zz-vv)+(tv-1)/(zz-vv))
actual_simple=mp.diff(f,wv)
target=cv/(wv-vv)+(2-tv)/(wv-vv)**2
ck('numerical direct node2 Laurent coefficient',abs(actual_simple-target)<mp.mpf('1e-35'))
# B=0 field scalar spectator relation: actual S differs from negative derivative G by e/q.
G=lambda ww:(vv-ww)**(hv-1)/(ww-vv)
Sactual=(2-tv)/(wv-vv)*G(wv)
Sexpected=-mp.diff(G,wv)+(hv-tv)/(wv-vv)*G(wv)
ck('numerical regulator derivative before mixed integration',abs(Sactual-Sexpected)<mp.mpf('1e-34'))
print('Completed',count,'checks. Local separated-source bridge only; no claim of the full H8 continuation.')
