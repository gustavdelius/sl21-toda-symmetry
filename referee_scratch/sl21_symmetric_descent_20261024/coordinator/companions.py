"""Fresh vector enumeration of all four quartet/bar families and source channels."""
import itertools,json
from pathlib import Path
import sympy as S
b,t,g,h,Q,R,H=S.symbols('b t g h Q R H',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I]);a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I]);T=a0+a2;d=-T/b**2;z=S.zeros(3,1)
dot=lambda p,q:S.factor(S.expand(p.dot(q)).subs(b**2,t))
vs={'A':(-d,z,0),'B':(d,z,1),'barA':(z,-d,0),'barB':(z,d,1),'0':(a0,a0,0),'2':(a2,a2,0),'V':(-T,-T,0)}
n=0;bad=[];records=[]
def ck(name,v):
 global n
 n+=1
 if not v:bad.append(name)
 print(('PASS ' if v else 'FAIL ')+name)
def degree(labels):
 power=S.Integer(0)
 for i,j in itertools.combinations(range(len(labels)),2):
  li,lj=labels[i],labels[j];a,c,_=vs[li];aa,cc,_=vs[lj]
  power-=dot(a,aa)+dot(c,cc)
  if li==lj and li in ('A','B'):power+=g-2/t
  if {li,lj}=={'0','2'}:power+=2*(h-t)
 derivative_count=sum(vs[k][2] for k in labels)
 nc=sum(k in ('A','B','barA','barB') for k in labels);ns=len(labels)-nc
 return S.factor(power-derivative_count+nc+2*ns-(1 if nc else 2))
for na_max,nb_max,bar in itertools.product([1,3],[1,3],['barA','barB']):
 if na_max==nb_max:continue
 for sources in [('0',),('2',),('V',),('0','0'),('2','2'),('0','2')]:
  for na,nb,nd in itertools.product(range(na_max+1),range(nb_max+1),range(2)):
   if na+nb+nd==0:continue
   labels=['A']*na+['B']*nb+[bar]*nd+list(sources);dg=degree(labels)
   old=(len(sources)==1 and ((len(labels)==2 and dg==0) or (len(labels)==3 and dg==-1)))
   initial=S.expand(S.cancel(dg.subs({g:Q+6,h:H+4,t:-2/(R+1)})))
   positive=all(c>=0 for c in S.Poly(initial,Q,R,H).coeffs()) and initial.subs({Q:0,R:0,H:0})>0
   ck('fixed proper face or positive initial union '+str((na_max,nb_max,bar,sources,na,nb,nd)),old or positive)
   if (na,nb,nd)==(na_max,nb_max,1):records.append({'quartet':f'A{na_max}B{nb_max}','bar':bar,'sources':sources,'full_degree':str(dg)})
# Core partial fraction exponents for a barred unary tail interacting with the other source.
for u,v,name in [(a0,a2,'barK0/Y2'),(a2,a0,'barK2/Y0')]:
 ck(name+' hol power h-1',S.simplify(-dot(u,v)+h-t-(h-1))==0)
 ck(name+' anti power h-2',S.simplify(-dot(d+u,v)+h-t-(h-2))==0)
# Endpoint test jet Taylor expansion with derivative contractions is automatic coefficient extraction.
r,rb=S.symbols('r rb');u0,u1,u2,f0,f1,f2=S.symbols('u0 u1 u2 f0 f1 f2')
N=u0+u1*r+u2*rb;test=f0+f1*r+f2*rb
ck('mixed first angular coefficient includes both test and field derivatives',S.expand(N*test).coeff(r,1).coeff(rb,1)==u1*f2+u2*f1)
Path(__file__).with_name('all_channel_degrees.json').write_text(json.dumps(records,indent=2)+'\n')
print('SUMMARY',n-len(bad),'PASS',len(bad),'FAIL')
if bad:raise SystemExit(1)
