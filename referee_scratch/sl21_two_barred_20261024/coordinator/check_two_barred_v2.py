"""Fresh primary-vector, full Wick, cycle and collision checks."""
import sympy as S
from itertools import combinations,product
from pathlib import Path
import json
out=Path(__file__).resolve().parent
b,t,gL,gR,h=S.symbols('b t gL gR h',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I])
a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I])
T=a0+a2;d=-T/b**2;zero=S.zeros(3,1)
dot=lambda v,w:S.factor(S.expand(v.dot(w)).subs(b**2,t))
checks=[]
def ck(name,value):
    good=bool(value);checks.append(good);print(('PASS ' if good else 'FAIL ')+name,flush=True)
vertices={'B':(d,zero,1),'A':(-d,zero,0),'barB':(zero,d,1),'barA':(zero,-d,0),'Y0':(a0,a0,0),'Y2':(a2,a2,0),'V':(-T,-T,0)}
def deg(labels):
    ncur=sum(x in ['A','B','barA','barB'] for x in labels)
    dimension=ncur+2*(len(labels)-ncur)-(1 if ncur else 2)
    exponent=S.Integer(0)
    for x,y in combinations(labels,2):
        l,r,_=vertices[x];ll,rr,_=vertices[y]
        exponent-=dot(l,ll)+dot(r,rr)
        if x==y and x in ['A','B']:exponent+=gL-2/t
        if x==y and x in ['barA','barB']:exponent+=gR-2/t
        if {x,y}=={'Y0','Y2'}:exponent+=2*(h-t)
    return S.factor(exponent+dimension-sum(vertices[x][2] for x in labels))
records=[]
for nl,nr in product(range(3),range(3)):
    if nl+nr==0:continue
    for src in [('Y0',),('Y2',),('Y0','Y0'),('Y2','Y2'),('Y0','Y2')]:
        actual=deg(['B']*nl+['barB']*nr+list(src))
        E=gL*nl*(nl-1)/2+gR*nr*(nr-1)/2
        n=nl+nr
        expected=E-n+1 if len(src)==1 else E-2*n+(2*h+1 if src==('Y0','Y2') else 5)
        ck('primary collision degree '+str((nl,nr,src)),S.simplify(actual-expected)==0)
        records.append({'left':nl,'right':nr,'sources':src,'degree':str(actual)})
ck('even-root four-A one-source degree',S.simplify(deg(['A','A','barA','barA','V'])-(gL+gR-3))==0)
ck('two-source full equal degree',deg(['B','B','barB','barB','Y2','Y2'])==gL+gR-3)
ck('two-source full unlike degree',deg(['B','B','barB','barB','Y0','Y2'])==gL+gR+2*h-7)
tstar=(7+S.sqrt(17))/4
omega=4/t+2*t-7
ck('unlike physical resonance candidate',S.simplify(omega.subs(t,tstar))==0)
ck('resonance lies strictly in prior physical patch',bool(2<tstar) and bool(tstar<3))
ck('unlike three-current unions positive at candidate',S.simplify((2/t+2*t-5).subs(t,tstar))>0)
# Exact scalar chiral H2 from actual derivative pair matching.
z1,z2,w0,w2,q,g=S.symbols('z1 z2 w0 w2 q g',nonzero=True)
r=z1-z2; f1=-1/r+q/(z1-w0)+q/(z1-w2);f2=1/r+q/(z2-w0)+q/(z2-w2)
H=f1*f2+1/r**2
L1=g/r-1/(z1-w0)-1/(z1-w2)
L2=-g/r-1/(z2-w0)-1/(z2-w2)
a=q*g-1
c=1-q-a*a+q*a
divform=q*q*(L1*L2+g/r**2)+q*a*((L1-L2)/r-2/r**2)+c/r**2
ck('complete scalar H2 divergence identity',S.factor(H-divform)==0)
cphysical=S.factor(c.subs(q,2/g-1))
ck('physical residual is nonzero, not pure external derivative',S.factor(cphysical+(1-g)*(2-g))==0)
ck('normalized edge shift cancels residual factors',S.factor(cphysical/((1-g)*(2-g))+1)==0)
# Cauchy determinant; no derivative numerator discarded in its stated scope.
det=1/((z1-w0)*(z2-w2))-1/((z1-w2)*(z2-w0))
ck('scalar Cauchy determinant sign',S.factor(det+(z1-z2)*(w0-w2)/((z1-w0)*(z1-w2)*(z2-w0)*(z2-w2)))==0)
# Two independent endpoint fusions: retain both ancestral current regulator shifts.
v=d+a2
power=-2*dot(v,v)+gL-2/t+gR-2/t
ck('double mixed endpoint base exponent',S.simplify(power-(gL+gR-2))==0)
ck('double endpoint base locally integrable physical patch',S.simplify(power.subs({gL:2/t,gR:2/t})+1-(4/t-1))==0)
ck('derivative endpoint degree agrees full primitive',S.simplify(power+1-2-(gL+gR-3))==0)
ck('endpoint degree lower bound at t=3',S.simplify((4/t-1).subs(t,3))==S.Rational(1,3))
# Loop exponents: a constant source-height gap provides two independent Laplace integrations.
G=S.symbols('G',real=True);k,span=S.symbols('k span',positive=True)
ck('source-height integrations lower primitive loop degree by two',S.simplify(k**(4-G)/(2*k)**2-k**(2-G)/4)==0)
ck('scaling of remaining Laplace integral',S.simplify((S.Symbol('u')/(2*span))**(2-G)/(2*span)/(S.Symbol('u')**(2-G)*(2*span)**(G-3)))==1)
(out/'collision_degrees.json').write_text(json.dumps(records,indent=2)+'\n')
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
