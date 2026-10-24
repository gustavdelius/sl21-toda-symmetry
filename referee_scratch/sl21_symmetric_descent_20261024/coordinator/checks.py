"""Fresh independent checks: full Wick matching, source masks, triangle and all union degrees."""
import sympy as S
import itertools,json
from pathlib import Path
out=Path(__file__).resolve().parent
checks=[]
def ck(name,yes):
 yes=bool(yes);checks.append(yes);print(('PASS ' if yes else 'FAIL ')+name,flush=True)
b,t,g,h,rho=S.symbols('b t g h rho',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I]);a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I]);d=-(a0+a2)/b**2;zero=S.zeros(3,1)
dot=lambda p,q:S.factor(S.expand(p.dot(q)).subs(b**2,t))
ck('primary Gram entries',dot(a0,a0)==-1 and dot(a2,a2)==-1 and dot(a0,a2)==1-t and dot(d,d)==-2/t)
# Generate derivative Wick matching recursively; every matching counted once.
zs=[S.Rational(-2),S.Rational(1,3),S.Rational(3)];a,w,kap=S.symbols('a w kap');cs=S.symbols('C0:3')
Fi=[cs[i]-sum(1/(zs[i]-zs[j]) for j in range(3) if j!=i)+1/(zs[i]-a)+kap/(zs[i]-w) for i in range(3)]
def match(ids):
 if not ids:return S.Integer(1)
 i,*rest=ids;answer=Fi[i]*match(tuple(rest))
 for j in rest:answer+=match(tuple(k for k in rest if k!=j))/(zs[i]-zs[j])**2
 return answer
actual=match((0,1,2));display=S.prod(Fi)+sum(Fi[3-i-j]/(zs[i]-zs[j])**2 for i,j in itertools.combinations(range(3),2))
ck('full three-C polynomial includes every derivative pairing',S.expand(actual-display)==0)
# Pole expansion includes arbitrary smooth coefficient/test jets, not just scalar residue.
q=S.symbols('q');f0,f1,f2,K,Kp,R,Rp=S.symbols('f0 f1 f2 K Kp R Rp')
physical=(t-1)*K/q**2-Kp/q
series=S.expand((f0+f1*q+f2*q*q/2)*(R+Rp*q)*physical)
ck('arbitrary moving-test double coefficient',S.expand(series.coeff(q,-2)-f0*R*(t-1)*K)==0)
ck('arbitrary moving-test simple coefficient',S.expand(series.coeff(q,-1)-(f0*(-R*Kp+(t-1)*Rp*K)+f1*(t-1)*R*K))==0)
ck('regulator derivative coefficient beta+1',S.expand(-R*Kp+(t-1)*Rp*K-(-(Rp*K+R*Kp)+t*Rp*K))==0)
# All one-source collision degrees computed from PRIMARY vector contractions.
records=[]
for ai in [a0,a2]:
 for na,nb,nd in itertools.product(range(2),range(4),range(2)):
  if na+nb+nd==0:continue
  vertices=[(-d,zero,'A')]*na+[(d,zero,'B')]*nb+[(zero,d,'D')]*nd+[(ai,ai,'Y')]
  power=S.Integer(0)
  for (l,r,L),(ll,rr,LL) in itertools.combinations(vertices,2):
   power-=dot(l,ll)+dot(r,rr)
   if L==LL=='B':power+=g-2/t
  degree=S.factor(power-nb-nd+(na+nb+nd+1))
  formula=g*nb*(nb-1)/2-(2/t)*na*nb+2*na-nb-nd+1
  ck('one-source union degree '+str((ai==a2,na,nb,nd)),S.simplify(degree-formula)==0)
  old=(na,nb,nd) in [(0,1,0),(0,0,1),(0,1,1)]
  if not old:
   # Choose g=G+6, rho=-(R+1), t=2/rho; exact strictly positive polynomial.
   GG,RR=S.symbols('GG RR',positive=True)
   pos=S.expand(formula.subs(t,-2/(RR+1)).subs(g,GG+6))
   ck('new one-source union vanishes in common domain '+str((ai==a2,na,nb,nd)),all(c>=0 for c in S.Poly(pos,GG,RR).coeffs()) and pos.subs({GG:0,RR:0})>0)
  if ai==a2:records.append({'A':na,'B':nb,'barB':nd,'degree':str(degree),'fixed_old_face':old})
(out/'one_source_faces.json').write_text(json.dumps(records,indent=2)+'\n')
# Mask derivatives: exact smooth radial surrogate; show all boundary terms by product rule.
eps,x,y=S.symbols('eps x y',positive=True);chi=S.Function('chi')
m1=chi((x*x+y*y)/eps**2);m2=chi(((x-1)**2+y*y)/eps**2)
ck('moving-mask first product derivative keeps both boundaries',S.simplify(S.diff(m1*m2,x)-(S.diff(m1,x)*m2+m1*S.diff(m2,x)))==0)
ck('moving-mask second derivative keeps mask intersection',S.simplify(S.diff(m1*m2,x,2)-(S.diff(m1,x,2)*m2+2*S.diff(m1,x)*S.diff(m2,x)+m1*S.diff(m2,x,2)))==0)
ck('hard-scale radial mask derivative has correct epsilon cost',S.simplify(S.diff(m1,x).subs({x:eps*S.Symbol('X'),y:eps*S.Symbol('Y')})*eps-2*S.Symbol('X')*S.Subs(S.Derivative(chi(S.Symbol('z')),S.Symbol('z')),S.Symbol('z'),S.Symbol('X')**2+S.Symbol('Y')**2))==0)
# Full unequal-source triangle, with the derivative-current terms retained.
p,q,cb=S.symbols('p q cb',nonzero=True);beta=t-1;rbar=q-p
original=(cb+1/p+beta/q)/(p*q)
reduced=(cb*(1/p-1/q)+1/p**2-beta/q**2)/rbar+(beta-1)*(1/p-1/q)/rbar**2
ck('unlike triangle reduction includes Cbar and both double poles',S.factor(original-reduced)==0)
# Homogeneities of extracted source pair multipliers and arbitrary finite differentiation.
N,HH=S.symbols('N HH',positive=True)
ck('source multipliers smooth enough after finite derivative loss',S.simplify((2*h-4-N).subs(h,(N+4)/2+HH))==2*HH)
# Mixed K0/K2 pair exponents and complete B3A degree.
K0=(d+a0,a0,'B');K2=(d+a2,a2,'B');B=(d,zero,'B');A=(-d,zero,'A')
ph=-dot(K0[0],K2[0])+g-2/t+h-t;pa=-dot(a0,a2)+h-t
ck('mixed-tail holomorphic exponent',S.simplify(ph-(g+h-3))==0)
ck('mixed-tail antiholomorphic exponent',S.simplify(pa-(h-1))==0)
ck('mixed-tail phase loses h dependence by even integer shift',S.simplify(ph-pa-(g-2))==0)
vs=[K0,K2,B,A];deg=S.Integer(0)
for i,j in itertools.combinations(range(4),2):
 l,r,L=vs[i];ll,rr,LL=vs[j];deg-=dot(l,ll)+dot(r,rr)
 if L==LL=='B':deg+=g-2/t
 if {i,j}=={0,1}:deg+=2*(h-t)
deg+=3-1
ck('full mixed-unary quartet positive degree',S.simplify(deg-(3*(g-2/t)+2*h-2))==0)
# Second-order coefficient identity for whole-word relation, no re-expansion per binary factor.
e0,e2=S.symbols('e0 e2');A00,A02,A22=S.symbols('A00 A02 A22')
poly=e0**2*A00+e0*e2*A02+e2**2*A22
ck('two labelled unary directions give cross coefficient exactly once',S.diff(poly,e0,e2).subs({e0:0,e2:0})==A02)
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
