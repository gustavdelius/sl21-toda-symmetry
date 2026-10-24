"""Independent even component Lax solve; Grassmann coefficients on the left."""
import sympy as s
from itertools import product
p=(0,0,1)
class G:
 def __init__(self,d=0):self.d=d if isinstance(d,dict) else {0:s.sympify(d)}
 def __add__(self,b):
  b=b if isinstance(b,G) else G(b)
  return G({k:s.expand(self.d.get(k,0)+b.d.get(k,0)) for k in self.d.keys()|b.d.keys()})
 __radd__=__add__
 def __neg__(self):return G({k:-v for k,v in self.d.items()})
 def __sub__(self,b):return self+-G(b) if not isinstance(b,G) else self+-b
 def __rsub__(self,b):return G(b)+-self
 def __mul__(self,b):
  b=b if isinstance(b,G) else G(b);out={}
  for i,u in self.d.items():
   for j,v in b.d.items():
    if i&j:continue
    sign=(-1)**sum(1 for a in range(4) for c in range(4) if a>c and (i>>a)&1 and (j>>c)&1)
    out[i|j]=out.get(i|j,0)+sign*u*v
  return G({k:s.expand(v) for k,v in out.items()})
 __rmul__=__mul__
 def __truediv__(self,b):return self*(1/b)
 def subs(self,d):return G({k:s.factor(v.subs(d)) for k,v in self.d.items()})
 def zero(self):return all(s.simplify(v)==0 for v in self.d.values())
 def __repr__(self):return str(self.d)
def zeros():return [[G() for j in range(3)] for i in range(3)]
def madd(a,b):return [[a[i][j]+b[i][j] for j in range(3)] for i in range(3)]
def neg(a):return [[-v for v in row] for row in a]
def mm(a,b):
 return [[sum((a[i][k]*b[k][j]*(-1)**((p[i]+p[k])*(p[k]+p[j])) for k in range(3)),G()) for j in range(3)] for i in range(3)]
def H(vec):
 u,v=vec;return [ (u+s.I*v)/s.sqrt(2),(-u+s.I*v)/s.sqrt(2),s.sqrt(2)*s.I*v]
rt=s.sqrt(2);roots=[s.Matrix([-1/rt,s.I/rt]),s.Matrix([rt,0]),s.Matrix([-1/rt,-s.I/rt])]
mu,z,A,B=s.symbols('mu z A B',nonzero=True)
a,b,c,d,e,f=s.symbols('a b c d e f')
Z=[G({1<<k:s.S.One}) for k in range(4)];pp=Z[:2];pm=Z[2:]
def dot(v,w):return sum((x*y for x,y in zip(v,w)),G())
S=[A,1/(A*B),B]
M=[[sum(mu*S[k]**2*roots[k][i]*roots[k][j] for k in (0,2)) for j in range(2)] for i in range(2)]
dmpp=[-dot(row,pm) for row in M];dppm=[dot(row,pp) for row in M]
Lp=zeros();Lm=zeros();DmLp=zeros();DpLm=zeros()
Lp[0][1]=G(a*mu/(A*B));Lp[1][0]=G(b*mu*z*A*B)
Lm[1][0]=G(d*mu/(A*B));Lm[0][1]=G(e*mu/z*A*B)
for k,i,j,mul in ((0,2,0,z),(2,1,2,1)):
 Lp[i][j]=c*s.sqrt(mu)*mul*S[k]*dot(roots[k],pp)
 DmLp[i][j]=c*s.sqrt(mu)*mul*S[k]*dot(roots[k],dmpp)
for k,i,j,mul in ((0,0,2,-1/z),(2,2,1,1)):
 Lm[i][j]=f*s.sqrt(mu)*mul*S[k]*dot(roots[k],pm)
 DpLm[i][j]=f*s.sqrt(mu)*mul*S[k]*dot(roots[k],dppm)
C=madd(madd(DpLm,neg(DmLp)),madd(mm(Lp,Lm),neg(mm(Lm,Lp))))
def main():
 eq=[]
 for i,j in product(range(3),repeat=2):
  if i!=j:
   for mask,v in C[i][j].d.items():
    if s.expand(v)!=0:
     print('OFF',i,j,mask,s.factor(v),flush=True);eq.append(s.factor(v))
 print('DIAGONAL',flush=True)
 for i in range(3):print(i,{k:s.factor(v) for k,v in C[i][i].d.items()},flush=True)
 print('SOLVE',s.solve(eq+[a*d+1,b*e+1,c*f+s.I],(a,b,c,d,e,f),dict=True),flush=True)

if __name__=="__main__":main()
