import sympy as s
from itertools import product
n=0

def chk(name,e):
 global n
 e=s.cancel(e)
 ok=e==0
 print(('PASS' if ok else 'FAIL')+' '+name,flush=True)
 if not ok:raise AssertionError((name,e))
 n+=1
x,y,r,q,g,C1,C2=s.symbols('x y r q g C1 C2')
Z=x-y
F1=C1-1/Z+1/x+q/(x-r);F2=C2+1/Z+1/y+q/(y-r)
H=F1*F2+1/Z**2
chk('complete unlike H2 scalar',H.subs({C1:0,C2:0})-q*((q-1)/((x-r)*(y-r))+1/(x*(y-r))+1/((x-r)*y)))
f1=C1-1/Z+q/x;f2=C2+1/Z+q/y;h=f1*f2+1/Z**2
chk('one-source full scalar',h.subs({C1:0,C2:0})-q*(q-1)/(x*y))
chk('one-source C1 coefficient',s.diff(h,C1).subs(C2,0)-(1/Z+q/y))
chk('one-source C2 coefficient',s.diff(h,C2).subs(C1,0)-(-1/Z+q/x))
D=1/(x*(y-r))-1/((x-r)*y)
chk('oriented equal-source Cauchy determinant',D-Z*r/(x*(x-r)*y*(y-r)))
# O1 O2 K/K from exact logarithmic derivatives, avoids symbolic fractional powers.
a=q*g-1;log1=g/Z-1/x;log2=-g/Z-1/y
O2K=C2-a/Z-q*log2
O1O2=(C1+a/Z)*O2K-q*(s.diff(O2K,x)+O2K*log1)
chk('full differential representative with residual',O1O2+(1-q)/Z**2-h)
D1,D2=s.symbols('D1 D2')
O2N=C2+q*D2-a/Z-q*(log2+D2)
O1O2N=(C1+q*D1+a/Z)*O2N-q*(s.diff(O2N,x)+O2N*(log1+D1))
chk('full independent normal first jets',O1O2N+(1-q)/Z**2-h)
k=s.symbols('k');A=1/r**2;B=s.I*k/r-1/r**2
P=(q-1)*A*B+B*B+A*A
chk('opposite-height residue polynomial',P-((s.I*k)**2/r**2+(q-3)*s.I*k/r**3+(3-q)/r**4))
chk('bulk current integral coefficient',g*(g-1)+(q-3)*g+3-q-(g-1)*(g+q-3))
t=s.symbols('t');c=(t-2)/2
chk('orthogonal-current independent coefficient',(q*(g-1)*(g+q-3)-(c-c*c*(g-1)*(g-4))).subs({q:t-1,g:2/t}))
ts=(7+s.sqrt(17))/4
assert s.simplify((4/t+2*t-7).subs(t,ts))==0
print('PASS physical full-union resonance');n+=1
cs=s.simplify(((t-1)*(2/t-1)*(2/t+t-4)).subs(t,ts))
assert cs!=0
print('PASS bulk residue coefficient nonzero '+str(cs));n+=1
for sources in [1,2]:
 for left,right in product(range(3),repeat=2):
  if left+right==0:continue
  E=g*(left*(left-1)+right*(right-1))/2
  omega=E-sources*(left+right)+(1 if sources==1 else 5)
  if sources==1 and left+right<=2:continue
  assert omega.subs(g,7)>0
  print('PASS high-g face '+str((sources,left,right))+' degree='+str(omega));n+=1
print('PASS total '+str(n))
