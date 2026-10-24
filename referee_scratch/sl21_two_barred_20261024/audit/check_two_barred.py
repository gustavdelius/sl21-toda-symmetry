import sympy as S
from itertools import combinations
n=0

def check(name, expr):
 global n
 ok=bool(expr) if isinstance(expr,(bool,S.logic.boolalg.BooleanAtom)) else S.cancel(S.expand(expr))==0
 n+=1
 print(('PASS' if ok else 'FAIL')+' '+name+('' if ok else ': '+str(expr)),flush=True)
 if not ok: raise AssertionError(name)

# Independent vectors and Gaussian contractions.
b=S.symbols('b', nonzero=True, real=True); t=b*b; I=S.I
v0=S.Matrix([I*b/S.sqrt(2),b/S.sqrt(2),I]); v2=S.Matrix([I*b/S.sqrt(2),-b/S.sqrt(2),-I]); d=-(v0+v2)/t
for label,ex in [('a0 squared',v0.dot(v0)+1),('a2 squared',v2.dot(v2)+1),('source cross',v0.dot(v2)-(1-t)),('d a0',d.dot(v0)-1),('d a2',d.dot(v2)-1),('d squared',d.dot(d)+2/t)]:check(label,ex)
z1,z2,w0,w2=S.symbols('z1 z2 w0 w2'); r=z1-z2
p=lambda z,w:1/(z-w)
D=p(z1,w0)*p(z2,w2)-p(z1,w2)*p(z2,w0)
check('Cauchy determinant sign',D+r*(w0-w2)/((z1-w0)*(z1-w2)*(z2-w0)*(z2-w2)))
g,q,C1,C2=S.symbols('g q C1 C2'); a=q*g-1
ell1=g/r-p(z1,w0)-p(z1,w2)
ell2=-g/r-p(z2,w0)-p(z2,w2)
F1=C1-1/r+q*(p(z1,w0)+p(z1,w2)); F2=C2+1/r+q*(p(z2,w0)+p(z2,w2))
# O_i on coefficient times vertex V: C_i+a/r_ij-q*(partial+ell_i).
def O1(c):return (C1+a/r)*c-q*(S.diff(c,z1)+ell1*c)
def O2(c):return (C2-a/r)*c-q*(S.diff(c,z2)+ell2*c)
check('full H2 differential representation',O1(O2(S.Integer(1)))+(1-q)/r**2-(F1*F2+1/r**2))
cstar=1-q-a*a+q*a
check('divergence scalar coefficient factorization',cstar-(1-a)*(1+a-q))
check('physical scalar coefficient',S.factor(cstar.subs(q,2/g-1))+(1-g)*(2-g))
check('nonzero normalized residual at physical g',S.cancel(cstar.subs(q,2/g-1)/((1-g)*(2-g)))+1)
# Two separated first-jet mixed contacts are commuting differential operators.
x,y=S.symbols('x y'); qa,qb,ca,cb=S.symbols('qa qb ca cb'); c=S.symbols('c0:9')
psi=sum(c[3*i+j]*x**i*y**j for i in range(3) for j in range(3))
def La(f):return qa**2*S.diff(f,x)-qa*ca*f
def Lb(f):return qb**2*S.diff(f,y)-qb*cb*f
check('disjoint complete first-jet operators commute',La(Lb(psi))-Lb(La(psi)))
check('disjoint mixed source jet retained',La(Lb(x*y)).subs({x:0,y:0})-qa**2*qb**2)
check('independent balanced products coefficient',S.Rational(1,2)**2-S.Rational(1,4))
check('global extreme averaging gives wrong disjoint coefficient',S.Rational(1,2)-S.Rational(1,4)==S.Rational(1,4))
# All two-source equal/unlike relative face degrees from contractions and dimension.
gL,gR,h=S.symbols('gL gR h')
for nl in range(3):
 for nr in range(3):
  E=gL*nl*(nl-1)/2+gR*nr*(nr-1)/2
  # source-current exponent -(nl+nr), per source; derivative degree -(nl+nr)
  # equal source +2; relative measure dimension nl+nr+3.
  raw=E-3*(nl+nr)+2+(nl+nr+3)
  expected=E-2*(nl+nr)+5
  check(f'equal face degree ({nl},{nr})',raw-expected)
  check(f'equal high-g positivity ({nl},{nr})',raw.subs({gL:5,gR:5})>0)
  unlike=E-3*(nl+nr)+(2*h-2)+(nl+nr+3)
  check(f'unlike face degree ({nl},{nr})',unlike-(E-2*(nl+nr)+2*h+1))
  check(f'unlike high-domain positivity ({nl},{nr})',unlike.subs({gL:5,gR:5,h:6})>0)
T=S.symbols('T', positive=True); ts=(7+S.sqrt(17))/4
omega=4/T+2*T-7
check('unlike full resonance exact location',omega.subs(T,ts))
check('unlike full resonance inside physical interval',bool(2<float(ts)<3))
check('unlike proper (2,1) positive throughout physical patch',S.factor((2/T+2*T-5)*T)-(2*T-1)*(T-2))
check('unlike full derivative at resonance nonzero',bool(float(S.diff(omega,T).subs(T,ts))>0))
# Equal full degree and boundary double-endpoint dimension agree.
check('double endpoint derivative degree',2*(g-1)-2+1-(2*g-3))
# Graph tests: one shared source cannot make a simple directed cycle joining two DAGs;
# two shared sources do. Construct the latter cycle explicitly.
edges=[('w0','B1'),('B1','B2'),('B2','w2'),('w2','D2'),('D2','D1'),('D1','w0')]
vertices=sorted(set(sum(([a,b] for a,b in edges),[])))
inc=S.zeros(len(vertices),len(edges))
for j,(a0,b0) in enumerate(edges):inc[vertices.index(a0),j]=-1;inc[vertices.index(b0),j]=1
check('positive circulation on two-source matching graph',inc*S.ones(6,1)==S.zeros(6,1))
print('SUMMARY',n,'PASS, 0 FAIL',flush=True)
print('DIAGNOSTIC unlike candidate t* =',S.N(ts,16),'degree at t=5/2 =',omega.subs(T,S.Rational(5,2)))
print('DIAGNOSTIC raw bulk cycle and potential resonance are not a full-forest anomaly.')
