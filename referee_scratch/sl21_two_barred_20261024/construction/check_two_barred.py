import sympy as s
from itertools import product
n=0

def chk(name,e):
 global n
 ok=s.simplify(e)==0
 print(('PASS' if ok else 'FAIL')+' '+name)
 if not ok:raise AssertionError((name,e))
 n+=1
z1,z2,w0,w2,q=s.symbols('z1 z2 w0 w2 q');C1,C2=s.symbols('C1 C2')
Z=z1-z2;p10=z1-w0;p12=z1-w2;p20=z2-w0;p22=z2-w2
F1=C1-1/Z+1/p10+q/p12;F2=C2+1/Z+1/p20+q/p22
H=s.expand(F1*F2+1/Z**2)
scalar=H.subs({C1:0,C2:0})
chk('complete unlike H2 scalar, including derivative pair',scalar-q*((q-1)/(p12*p22)+1/(p10*p22)+1/(p12*p20)))
# One source arbitrary coupling: all three Wick matchings survive.
w=s.symbols('w');f1=C1-1/Z+q/(z1-w);f2=C2+1/Z+q/(z2-w)
h=s.expand(f1*f2+1/Z**2)
chk('one-source H2 scalar',(h.subs({C1:0,C2:0})-q*(q-1)/((z1-w)*(z2-w))))
chk('one-source H2 coefficient C1',h.coeff(C1,1).coeff(C2,0)-(1/Z+q/(z2-w)))
chk('one-source H2 coefficient C2',h.coeff(C2,1).coeff(C1,0)-(-1/Z+q/(z1-w)))
# Exact equal-source Cauchy determinant, including orientation sign.
D=1/(p10*p22)-1/(p12*p20)
chk('equal-source determinant',D+(Z*(w0-w2))/(p10*p12*p20*p22))
# Full differential representative, keep residual (1-q)/Z^2.
g,D1,D2=s.symbols('g D1 D2')
K=Z**g/((z1-w)*(z2-w))
a=q*g-1
O2K=(C2-a/Z)*K-q*s.diff(K,z2)
O1O2=(C1+a/Z)*O2K-q*s.diff(O2K,z1)
chk('one-source full differential identity',(O1O2+(1-q)*K/Z**2-K*h)/K)
# Full field jets: derivative normal exponential is canceled by +q Di.
N=s.exp(D1*z1+D2*z2)
def O1(A):return (C1+q*D1+a/Z)*A-q*s.diff(A,z1)
def O2(A):return (C2+q*D2-a/Z)*A-q*s.diff(A,z2)
chk('same identity retains independent normal first jets',(O1(O2(K*N))+(1-q)*K*N/Z**2-K*N*h)/(K*N))
# Two-source hol current integral: exact residue polynomial before k integration.
k,r=s.symbols('k r');A=1/r**2;B=s.I*k/r-1/r**2
P=s.expand((q-1)*A*B+B*B+A*A)
chk('opposite-height residue polynomial',P-((s.I*k)**2/r**2+(q-3)*s.I*k/r**3+(3-q)/r**4))
g=s.symbols('g')
coef=g*(g-1)+(q-3)*g+3-q
chk('full unlike current-integral coefficient',coef-(g-1)*(g+q-3))
t=s.symbols('t',positive=True);c=(t-2)/2
chk('orthogonal-current independent coefficient comparison',(q*(g-1)*(g+q-3)-(c-c*c*(g-1)*(g-4))).subs({q:t-1,g:2/t}))
tstar=(7+s.sqrt(17))/4
chk('new full-union physical resonance',(4/t+2*t-7).subs(t,tstar))
cs=s.simplify(((t-1)*(2/t-1)*(2/t+t-4)).subs(t,tstar))
assert cs!=0
print('PASS unlike bulk coefficient nonzero at candidate resonance: '+str(cs));n+=1
# Proper source-containing face degrees, symmetric current counts <=2 each.
for sources in [1,2]:
 for left,right in product(range(3),repeat=2):
  if left+right==0:continue
  E=g*(left*(left-1)+right*(right-1))/2
  omega=E-sources*(left+right)+(1 if sources==1 else 5)
  if sources==1 and left+right<=2:continue
  assert omega.subs(g,7)>0
  print('PASS high-g face '+str((sources,left,right))+' degree='+str(omega));n+=1
print('PASS total '+str(n))
