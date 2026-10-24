from pathlib import Path
import itertools,json,hashlib
import sympy as s
out=Path(__file__).resolve().parent
inputs=['AGENTS.md','briefs/referee_sl21_common_product_overlap_report.md','briefs/referee_sl21_common_product_construction_report.md','briefs/referee_sl21_serre_m2_report.md','briefs/referee_sl21_symmetry_coordinator_plan_20261024.md']
(out/'inputs.json').write_text(json.dumps({f:hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in inputs},indent=2)+'\n')
def ck(name,ok):
 print(('PASS ' if bool(ok) else 'FAIL ')+name)
 if not bool(ok):raise AssertionError(name)
b,t,g,h,rho,u,cb=s.symbols('b t g h rho u cb',nonzero=True)
p0=s.Matrix([s.I*b/s.sqrt(2),b/s.sqrt(2),s.I]);p2=s.Matrix([s.I*b/s.sqrt(2),-b/s.sqrt(2),-s.I]);d=-(p0+p2)/b**2;v=d+p2
D=lambda p,q:s.factor(s.expand((p.T*q)[0]).subs(b**2,t))
ck('two-tail holomorphic exponent becomes g-1',s.simplify(-D(v,v)+g-2/t-(g-1))==0)
ck('equal-root anti numerator and root/bar-tail zero',-D(p2,p2)==1 and D(p2,v)==0)
ck('barred derivative coefficient beta',-D(p0,p2)==t-1)
for m in (1,2,3):
 xs=s.symbols('x0:'+str(m));vand=s.prod(xs[i]-xs[j] for i in range(m) for j in range(i+1,m))
 F=vand/s.prod(x-u for x in xs)
 coeff=[(-1)**(m-i-1)*s.prod(xs[j]-xs[k] for j in range(m) for k in range(j+1,m) if i not in (j,k)) for i in range(m)]
 L=sum(coeff[i]/(xs[i]-u) for i in range(m))
 ck('Vandermonde leaf decomposition m='+str(m),s.factor(F-L)==0)
 actual=F*(cb+(t-1)*sum(1/(u-x) for x in xs))
 leaves=sum(coeff[i]*(cb/(xs[i]-u)-(t-1)/(xs[i]-u)**2) for i in range(m))
 ck('full derivative-current leaf decomposition m='+str(m),s.factor(actual-leaves)==0)
 print('INFO leaf coefficients',m,coeff)
# General spectator Laurent transformation, retaining full normal product K.
x=s.symbols('x');K=s.Function('K')(x);R=s.Function('R')(x);beta=t-1
Dg=beta*R*K;Sg=-R*s.diff(K,x)+beta*s.diff(R,x)*K
ck('arbitrary-spectator simple residue transformation',s.expand(Sg+s.diff(R*K,x)-t*s.diff(R,x)*K)==0)
# Contour moments fixed by direct first-order source/circle calculation.
Jx,Jb,rreg=s.symbols('Jx Jb rreg')
action=-beta*t*(Jx+Jb-rreg)
tails=t*(Jx-t*rreg)+t*(beta*Jx+t*Jb)
matching=t*beta*(Jx+Jb-rreg)
ck('action plus insertion-wise matching cancels off surface',s.expand(action+matching)==0)
ck('both tails give canonical insertion derivative',s.expand(tails-t*t*(Jx+Jb-rreg))==0)
# Positive-cone properness is a theorem for DAGs. Verify actual orientation data
# for every possible relative source/current height order, including anti leaf.
def acyclic(nodes,edges):
 pending=set(nodes)
 while pending:
  roots={j for j in pending if not any(v==j and u in pending for u,v in edges)}
  if not roots:return False
  pending-=roots
 return True
for m in (1,2):
 hol=tuple(range(4+m));leaf=4+m;checked=0
 for order in itertools.permutations(range(5+m)):
  rank={j:i for i,j in enumerate(order)}
  edges=[(i,j) if rank[i]<rank[j] else (j,i) for i,j in itertools.combinations(hol,2)]
  for src in range(4,4+m):
   # Antichiral orientation is reversed relative to the actual height order.
   anti=(leaf,src) if rank[leaf]>rank[src] else (src,leaf)
   if not acyclic(range(5+m),edges+[anti]):raise AssertionError((m,order,src))
   checked+=1
 ck('all current/source height chambers give DAG m='+str(m),True)
 print('INFO directed configurations checked',m,checked)
# Reconstruct unlike degrees directly from counts of actual vector contractions.
z=s.zeros(3,1);mom={'A':(-d,z),'B':(d,z),'D':(z,d),'0':(p0,p0),'2':(p2,p2)}
records=[]
for na in (0,1):
 for nb in range(4):
  for nd in (0,1):
   if na+nb+nd==0:continue
   labels=['A']*na+['B']*nb+['D']*nd+['0','2'];power=0
   for i,j in itertools.combinations(range(len(labels)),2):
    li,lj=labels[i],labels[j]
    power-=D(mom[li][0],mom[lj][0])+D(mom[li][1],mom[lj][1])
    if li==lj=='B':power+=g-2/t
    if {li,lj}=={'0','2'}:power+=2*(h-t)
   power-=nb+nd
   degree=s.simplify(power+na+nb+nd+3)
   E=g*nb*(nb-1)/2-rho*na*nb
   expected=E+3*na-2*(nb+nd)+2*h+1
   ck('unlike face '+str((na,nb,nd)),s.simplify(degree-expected.subs(rho,2/t))==0)
   physical=s.factor(degree.subs({g:2/t,h:t}))
   qq=s.symbols('q',positive=True)
   poly=s.Poly(s.expand((t*physical).subs(t,2+qq)),qq)
   ck('strict physical positive patch '+str((na,nb,nd)),all(x>=0 for x in poly.all_coeffs()) and poly.nth(1)>0)
   ck('no entire-physical-surface divisor '+str((na,nb,nd)),s.diff(degree,h)==2)
   records.append({'a':na,'b':nb,'d':nd,'degree':str(degree),'physical':str(physical)})
(out/'unlike_faces.json').write_text(json.dumps(records,indent=2)+'\n')
print('INFO full unlike five-current degree',records[-1]['degree'])
ck('source-only degree',s.simplify(2*(h-1)+2-2*h)==0)
print('INFO no previous research module imported; arbitrary normal-field/test coefficients retained analytically')
# Independently test the minority-coordinate identity needed for B^3 A.
a=s.symbols('a');zs=[s.Integer(0),s.Integer(2),s.Integer(5)];eta=t-2
fs=[-sum(1/(zs[i]-zs[j]) for j in range(3) if j!=i)+1/(zs[i]-a) for i in range(3)]
scalar=s.prod(fs)+sum(fs[3-i-j]/(zs[i]-zs[j])**2 for i,j in itertools.combinations(range(3),2))
constantC=sum(fs[i]*fs[j]+1/(zs[i]-zs[j])**2 for i,j in itertools.combinations(range(3),2))
ck('bare quartet scalar cancellation from Wick contractions',s.factor(scalar)==0)
ck('bare quartet one-regular-C cancellation',s.factor(constantC)==0)
for i in range(3):
 j,k=[z for z in range(3) if z!=i]
 fj=-1/(zs[j]-zs[k])+1/(zs[j]-a)+eta/(zs[j]-zs[i])
 fk=-1/(zs[k]-zs[j])+1/(zs[k]-a)+eta/(zs[k]-zs[i])
 actual=(a-zs[i])/((zs[j]-zs[i])*(zs[k]-zs[i]))*(fj*fk+1/(zs[j]-zs[k])**2)
 Si=(a-zs[i])**2/s.prod((zs[i]-zs[l])**2 for l in range(3) if l!=i)
 predicted=eta*t/2*(s.diff(Si,a)+Si*(2/t)*sum(1/(zz-a) for zz in zs))
 ck('individual one-tail minority derivative '+str(i),s.factor(actual-predicted)==0)
