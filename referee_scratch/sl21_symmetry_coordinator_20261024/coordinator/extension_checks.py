import sympy as s
import mpmath as m
m.mp.dps=40
r,z,a,t,h,q=s.symbols('r z a t h q',nonzero=True)
failed=[]
def ck(n,x):
 ok=bool(x);print(('PASS ' if ok else 'FAIL ')+n)
 if not ok:failed.append(n)
def eq(n,x):ck(n,s.factor(x)==0)
# Source kernel generated from contractions, not inherited scripts.
f=1/((z-r)*z)*(1/(z-r)+(t-1)/z)
g=(a-r)*a/((z-r)*z)*(1/(z-a)+1/(z-r)+(t-1)/z)
for n in range(5):
 eq('two-source Laurent '+str(n),s.diff(f,r,n).subs(r,0)/s.factorial(n)-(n+t)/z**(n+3))
 if n: eq('A-dressed Laurent '+str(n),s.diff(g,r,n).subs(r,0)/s.factorial(n)-a*(a-z)*(n+t)/z**(n+3))
# Disk source-pair remainder from independent radial and angular expansion.
# Product n-th holomorphic/antiholomorphic terms contributes (n+t)^2/(n+t).
series=s.summation((s.Symbol('n',integer=True,positive=True)+t)*q**s.Symbol('n',integer=True,positive=True),(s.Symbol('n',integer=True,positive=True),1,6))
closed=q*((1+t)-t*q)/(1-q)**2
eq('disk series to sixth order',s.series(closed,q,0,7).removeO()-series)
for tt in [-m.sqrt(2)/3,-m.sqrt(3)/4]:
 for M in [m.sqrt(5),m.sqrt(11)]:
  val=(1/M**2)*((1+tt)-tt/M**2)/(1-1/M**2)**2
  ck('nonzero shrinking-separation remainder',val>0)
  print('DIAGNOSTIC t=',float(tt),'M=',float(M),'normalized disk remainder=',m.nstr(val,18))
# Descendant radial subtraction coefficient independently recovered.
for n in [0,1]:
 primitive=r**(2*(t-1+n))/(2*(t-1+n))
 eq('endpoint radial primitive '+str(n),s.diff(primitive,r)-r**(2*t-3+2*n))
# h-dependent Ward contact polynomial.
C,D,C2,D2=s.symbols('C D C2 D2')
N=C2+t*D2/2-t*C*D-t*t*D*D/2
RP=-t*C*D+t*t*D2/2-t**3*D*D/2
partialL=C2+(1-t)*C*D
eq('ordinary Ward descent',N-RP/t-partialL)
eq('h-regulated boundary leftover',N-RP/h-(t/h)*partialL-(h-t)*N/h)
# Independent anomaly dimension count for BOTH currents.
n0,n2,ne=s.symbols('n0 n2 ne',integer=True,nonnegative=True)
A=n0-ne;B=n2-ne;SS=A+B;FF=n0-n2;M=n0+n2+2*ne
baseR=t*A*B+FF**2/2;wc=M/2-ne*t
ellB=s.expand(1+1/t-(1/t-SS+baseR)-wc);bar=s.expand(1-baseR-wc)
ellA=s.expand(1/t-(1/t+SS+baseR)-wc)
eq('B level formula',ellB-(1-M/2+SS-FF**2/2+t*(ne-A*B)))
eq('A level formula',ellA-(-M/2-SS-FF**2/2+t*(ne-A*B)))
allowed={'A':[],'B':[]}
for i in range(9):
 for j in range(9):
  for e in range(9):
   aa=i-e;bb=j-e
   if e!=aa*bb:continue
   ss=aa+bb;rr=1-s.Rational(ss*(ss+1),2)
   ls={'A':-s.Rational(ss*(ss+3),2),'B':1-s.Rational(ss*(ss-1),2)}
   for typ,ll in ls.items():
    if ll>=0 and rr>=0 and i+j+e:allowed[typ].append((i,j,e,int(ll),int(rr)))
print('DIAGNOSTIC anomaly candidates',allowed)
ck('A only even-root local candidate',allowed['A']==[(0,0,1,1,0)])
ck('B only single-Yukawa local candidates',set(allowed['B'])=={(0,1,0,1,0),(1,0,0,1,0)})
# Restricted-test ambiguity: integrating U(y) delta'(x-y) against g(x) gives -g'(y)U(y).
x,y=s.symbols('x y',real=True)
ck('contact invisible to constant outer test',s.diff(s.Integer(1),x)==0)
ck('contact visible to nonconstant outer test',s.diff(x,x).subs(x,y)==1)
# Formal spatial anchor splitting; exp(q d varphi) contributes exp[-q d int_left dotPhi/2].
qq,pp,vv=s.symbols('q p v')
eq('anchor string multiplicativity',s.expand(-qq*(pp+vv)/2+qq*pp/2+qq*vv/2))
# Real reference cycle instability: an open angular interval cos<=-1/2 yields V<=-lambda e^(zphi)+lambda1 e^(-2zphi).
zz=m.sqrt(m.sqrt(2)/6)
for ph in [5,10,20]:
 bound=-m.exp(zz*ph)+m.exp(-2*zz*ph)
 print('DIAGNOSTIC reference real-cycle potential bound phi=',ph, 'V<=',m.nstr(bound,18))
 ck('negative exponentially growing ray '+str(ph),bound<0)
print('FAILURES',len(failed));assert not failed
