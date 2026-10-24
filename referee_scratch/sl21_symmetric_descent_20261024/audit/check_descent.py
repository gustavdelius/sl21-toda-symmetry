import sympy as S
from itertools import permutations,combinations
from pathlib import Path
import hashlib,json
out=Path(__file__).resolve().parent;root=out.parents[2]
inputs=['AGENTS.md','briefs/referee_sl21_symmetry_coordinator_plan_20261024.md','briefs/referee_sl21_interaction_map_report.md','briefs/referee_sl21_interaction_map_audit_report.md']
(out/'input_hashes.json').write_text(json.dumps({x:hashlib.sha256((root/x).read_bytes()).hexdigest() for x in inputs},indent=2))
count=0
def check(name,expr):
 global count
 expr=S.factor(S.cancel(expr));ok=expr==0
 print(('PASS' if ok else 'FAIL')+' '+name+' residual='+str(expr));count+=1
 if not ok:raise AssertionError(name)
p,q,C,beta=S.symbols('p q C beta');rbar=q-p
base=(C+1/p+beta/q)/(p*q)
resolved=(C*(1/p-1/q)+1/p**2-beta/q**2)/rbar+(beta-1)*(1/p-1/q)/rbar**2
check('unlike triangle full current polynomial',base-resolved)
k=S.symbols('k')
check('equal triangle full current polynomial',rbar/(p*q)*(C+k/p+k/q)-(C*(1/p-1/q)+k*(1/p**2-1/q**2)))
b=S.symbols('b',positive=True);t=b*b;rho=2/t;g,h=S.symbols('g h');I=S.I
a0=S.Matrix([I*b/S.sqrt(2),b/S.sqrt(2),I]);a2=S.Matrix([I*b/S.sqrt(2),-b/S.sqrt(2),-I]);d=-(a0+a2)/t;zero=S.zeros(3,1)
for n,e in [('a0sq',a0.dot(a0)+1),('a2sq',a2.dot(a2)+1),('d.a0',d.dot(a0)-1),('d.a2',d.dot(a2)-1),('a0.a2',a0.dot(a2)-(1-t)),('bar coefficients node0',-a0.dot(a0)-1),('bar coefficients node2',-a0.dot(a2)-(t-1))]:check(n,e)
K0=d+a0;K2=d+a2
check('mixed tail hol exponent',-K0.dot(K2)+(g-rho)+(h-t)-(g+h-3))
check('mixed tail anti exponent',-a0.dot(a2)+(h-t)-(h-1))
check('mixed tail total exchange exponent h cancels',(g+h-3)-(h-1)-(g-2))
check('tail source unlike hol exponent',-K2.dot(a0)+(h-t)-(h-2))
check('tail source unlike anti exponent',-a2.dot(a0)+(h-t)-(h-1))
check('tail source equal anti polynomial',-a2.dot(a2)-1)
check('tail source equal hol nonsingular',-K2.dot(a2))
# Physical vector dot products reconstructed directly; no old degree code imported.
def rec(H,A,cur=None,der=0,src=None,dim=1):return dict(H=H,A=A,cur=cur,der=der,src=src,dim=dim)
B=rec(d,zero,'B',1);A=rec(-d,zero,'A');D=rec(zero,d,'D',1)
Y0=rec(a0,a0,None,0,0,2);Y2=rec(a2,a2,None,0,2,2)
T0=rec(K0,a0,'B',0,0);T2=rec(K2,a2,'B',0,2)
def degree(vertices):
 e=sum(x['dim'] for x in vertices)-1-sum(x['der'] for x in vertices)
 for x,y in combinations(vertices,2):
  e-=x['H'].dot(y['H'])+x['A'].dot(y['A'])
  if x['cur']==y['cur'] and x['cur'] in ('A','B'):e+=g-rho
  if {x['src'],y['src']}=={0,2}:e+=2*(h-t)
 return S.factor(e)
for species,amax,bmax in [('B',1,3),('A',3,1)]:
 for aa in range(amax+1):
  for bb in range(bmax+1):
   for dd in range(2):
    if aa+bb+dd==0:continue
    vv=[A]*aa+[B]*bb+[D]*dd+[Y0,Y2]
    E=g*(aa*(aa-1)+bb*(bb-1))/2-rho*aa*bb
    formula=E+3*aa-2*(bb+dd)+2*h+1
    check('new unlike degree '+species+str((aa,bb,dd)),degree(vv)-formula)
    physical=S.factor(degree(vv).subs({g:rho,h:t}))
    x=S.symbols('x',nonnegative=True)
    poly=S.Poly(S.expand((t*physical).subs(b*b,2+x)),x)
    assert all(c>=0 for c in poly.all_coeffs()) and any(c>0 for c in poly.all_coeffs())
    print('PASS positive physical unlike face '+species+str((aa,bb,dd))+' t*degree='+str(poly.as_expr()));count+=1
# Mixed double-unary proper and full faces.
for name,vv,want in [('pair',[T0,T2],g+2*h-3),('B triple',[T0,T2,B],3*g+2*h-5),('A triple',[T0,T2,A],g-2*rho+2*h),('full quartet',[T0,T2,B,A],3*(g-rho)+2*h-2)]:check('mixed double unary '+name,degree(vv)-want)
check('equal sources whole five degree',degree([B,B,B,A,D,Y2,Y2])-3*(g-rho))
# Source-to-tail weights retain full physical derivative-current polynomials.
z,v,a,c=S.symbols('z v a c');eta=t-2
fz=-1/(z-v)+1/(z-a)+eta/(z-c);fv=1/(z-v)+1/(v-a)+eta/(v-c)
Sp=(a-c)**2/((c-z)**2*(c-v)**2)
logRa=rho*(1/(z-a)+1/(v-a)+1/(c-a))
check('K2 marginal scalar primitive',(a-c)/((z-c)*(v-c))*(fz*fv+1/(z-v)**2)-(eta/rho)*(S.diff(Sp,a)+Sp*logRa))
f0=-1/(z-v)+1/(z-a);f1=1/(z-v)+1/(v-a)
check('K0 scalar zero',f0*f1+1/(z-v)**2)
a1,a2p,a3=S.symbols('a1 a2p a3');SA=(z-a1)**2/((a2p-a1)**2*(a3-a1)**2)
logRz=-rho*(1/(z-a1)+1/(z-a2p)+1/(z-a3))
check('FA marginal scalar primitive',SA*(1/(z-a1)+1/(z-a2p)+1/(z-a3)-t/(z-a1))+(S.diff(SA,z)+SA*logRz)/rho)
# Quantum shuffle with ALL majority slots decorated distinctly; proves polarized
# coefficient identities for bare, unary, double-unary assignments together.
s=S.symbols('s',nonzero=True);cs=1+s+1/s
letters=['X0','X2','B','A'];major=['X0','X2','B'];coef=[1,-cs,cs,-1]
def weight(word,target):
 pos={v:i for i,v in enumerate(target)};ans=1
 for i in range(4):
  for j in range(i+1,4):
   if pos[word[i]]>pos[word[j]]:ans*=s if (word[i]=='A')==(word[j]=='A') else 1/s
 return ans
for target in permutations(letters):
 total=0
 for kA in range(4):
  for maj in permutations(major):
   word=list(maj);word.insert(3-kA,'A');total+=coef[kA]*weight(word,target)
 check('fully decorated shuffle '+''.join(target),total)
# Forest factorials: one source assigned once, no same-current double source.
lam0,lam2,f0,f2=S.symbols('lam0 lam2 f0 f2')
exp2=S.expand((lam0*f0+lam2*f2)**2/2)
check('distinct source action factor',exp2.coeff(lam0).coeff(lam2)-f0*f2)
check('equal source action factor',exp2.coeff(lam2,2)-f2**2/2)
print('TOTAL',count,'PASS')
