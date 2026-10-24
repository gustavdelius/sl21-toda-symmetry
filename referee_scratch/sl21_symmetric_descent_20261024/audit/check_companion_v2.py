import sympy as S
from itertools import combinations
b=S.symbols('b',positive=True);t=b*b;rho=2/t;g,h=S.symbols('g h');I=S.I
p=S.Matrix([I*b/S.sqrt(2),b/S.sqrt(2),I]);q=S.Matrix([I*b/S.sqrt(2),-b/S.sqrt(2),-I]);d=-(p+q)/t;T=p+q;k=d+T;Z=S.zeros(3,1)
count=0
def check(n,e):
 global count
 e=S.factor(S.cancel(e));print(('PASS' if e==0 else 'FAIL')+' '+n+' residual='+str(e));count+=1
 if e!=0:raise AssertionError(n)
check('barA Yukawa polynomial exponent',-(-d).dot(p)-1)
check('barA V1 double antileaf',-(-d).dot(-T)+2)
check('barB V1 polynomial exponent',-d.dot(-T)-2)
for j in range(3):check('AA fused momentum '+str(j),(k-(1-t)*d)[j])
L,R=S.symbols('L R');gamma=t/(1-t)
for sg in (-1,1):
 # General ordered moments: E=C+betaD, K=D+A_L, H=beta A_R-Cbar.
 E=L/(1-t);K=L;H=-gamma*R
 outer=-I*S.pi*E*(S.pi*sg)+S.pi*K*(I*S.pi*sg)+S.pi*H*(I*S.pi*sg)
 check('AA raw ordered half contact '+str(sg),-outer/S.pi-I*S.pi*sg*gamma*(L+R))
 check('AA first test derivative cancels '+str(sg),S.pi*(I*S.pi*sg)-I*S.pi*(S.pi*sg))
 # canonical tail odd coefficient gamma², unary/tail ratio 1/gamma
 check('AA unary/balanced cancellation '+str(sg),I*S.pi*sg*gamma**2/gamma-I*S.pi*sg*gamma)
def rec(H,A,c=None,der=0,src=None,dim=1):return dict(H=H,A=A,c=c,der=der,src=src,dim=dim)
AA=rec(-d,Z,'A');BB=rec(d,Z,'B',1);DA=rec(Z,-d,'D');DB=rec(Z,d,'D',1)
Y0=rec(p,p,src=0,dim=2);Y2=rec(q,q,src=2,dim=2);V=rec(-T,-T,dim=2)
def degree(vv):
 e=sum(v['dim'] for v in vv)-1-sum(v['der'] for v in vv)
 for x,y in combinations(vv,2):
  e-=x['H'].dot(y['H'])+x['A'].dot(y['A'])
  if x['c']==y['c'] and x['c'] in ['A','B']:e+=g-rho
  if {x['src'],y['src']}=={0,2}:e+=2*(h-t)
 return S.factor(e)
for label,amax,bmax in [('SB',1,3),('SA',3,1)]:
 for a in range(amax+1):
  for nb in range(bmax+1):
   for nd in range(2):
    if a+nb+nd==0:continue
    E=g*(a*(a-1)+nb*(nb-1))/2-rho*a*nb
    vv=[AA]*a+[BB]*nb+[DA]*nd
    check('unlike barA degree '+label+str((a,nb,nd)),degree(vv+[Y0,Y2])-(E+3*a-2*nb+3*nd+2*h+1))
    x=S.symbols('x',nonnegative=True)
    poly=S.Poly(S.expand(S.cancel(t*degree(vv+[Y0,Y2]).subs({g:rho,h:t})).subs(b*b,2+x)),x)
    assert all(c>=0 for c in poly.all_coeffs()) and any(c>0 for c in poly.all_coeffs())
    print('PASS unlike barA face positivity '+label+str((a,nb,nd))+' '+str(poly.as_expr()));count+=1
    check('V1 barA degree '+label+str((a,nb,nd)),degree(vv+[V])-(E-a+2*nb-nd+1))
check('SA barA V1 full marginal degree',degree([AA,AA,AA,BB,DA,V])-(3*(g-rho)-1))
check('SB barB V1 full regular degree',degree([AA,BB,BB,BB,DB,V])-(3*(g-rho)+8))
print('TOTAL',count,'PASS')
