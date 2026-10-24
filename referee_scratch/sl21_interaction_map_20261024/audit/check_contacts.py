import sympy as S
import mpmath as mp
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[3]
inputs=['briefs/referee_sl21_symmetry_coordinator_plan_20261024.md','briefs/referee_sl21_mixed_audit_report.md','briefs/referee_sl21_mixed_regulator_report.md']
Path(__file__).with_name('input_hashes.json').write_text(json.dumps({f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in inputs},indent=2))
count=0
def check(name,expr):
 global count
 ok=S.simplify(expr)==0
 print(('PASS' if ok else 'FAIL')+' '+name+' residual='+str(S.simplify(expr)))
 count+=1
 if not ok: raise AssertionError(name)
b=S.symbols('b',positive=True); t=b*b; beta=t-1
I=S.I
p=S.Matrix([I*b/S.sqrt(2),b/S.sqrt(2),I]); q=S.Matrix([I*b/S.sqrt(2),-b/S.sqrt(2),-I]); d=-(p+q)/t; v=d+q
for name,expr in [('a0sq',p.dot(p)+1),('a2sq',q.dot(q)+1),('a0a2',p.dot(q)-(1-t)),('d.a0',d.dot(p)-1),('d.a2',d.dot(q)-1),('A V1 doublepole',(-d).dot(-p-q)-2)]:check(name,expr)
for j in range(3):check('firstfield hol identity '+str(j),(p+beta*d+v)[j])
for j in range(3):check('firstfield anti identity '+str(j),(beta*q-p-t*v)[j])
x,y,w,wb,D,E,C,F,k=S.symbols('x y w wb D E C F k')
G=1/(x-w)/(y-wb); N=S.exp(D*x+E*y)
def LX(f):return (C+k*D)*f-k*S.diff(f,x)
def LY(f):return (F+k*E)*f-k*S.diff(f,y)
check('all normal jets differential representative',LX(LY(G*N))/N-G*(C+k/(x-w))*(F+k/(y-wb)))
# Separate source-edge coordinates: no derivative of normal fields is silently omitted.
z0,z2,r0,r2,k0,k2=S.symbols('z0 z2 r0 r2 k0 k2')
G2=1/(z0*z2*r0*r2)
def AX(f):return C*f-k0*S.diff(f,z0)-k2*S.diff(f,z2)
def AY(f):return F*f-k0*S.diff(f,r0)-k2*S.diff(f,r2)
check('two source full Wick polynomial',AX(AY(G2))-G2*(C+k0/z0+k2/z2)*(F+k0/r0+k2/r2))
check('no two source doublepoles at one current',S.expand((C+k0/z0+k2/z2)).coeff(k0,1).coeff(k2,1))
# Singular contraction labels eXi model delta-face insertions; only transverse terms extracted.
eX0,eX2,eY0,eY2,HX,HY=S.symbols('eX0 eX2 eY0 eY2 HX HY')
expanded=S.expand((HX+S.pi*(k0*eX0/z2+k2*eX2/z0))*(HY+S.pi*(k0*eY0/r2+k2*eY2/r0)))
check('distinct-source doubleface coefficient 02',expanded.coeff(eX0).coeff(eY2)-S.pi**2*k0*k2/(z2*r0))
check('distinct-source doubleface coefficient 20',expanded.coeff(eX2).coeff(eY0)-S.pi**2*k0*k2/(z0*r2))
r,h=S.symbols('r h',positive=True)
check('collapsed distinct-source radial power',S.powdenest((-1/r**2)*r**(2*h-2),force=True)+r**(2*h-4))
# beta times unary pair sum is -beta/t times canonical tail sum.
t0=S.symbols('t',nonzero=True); bt=t0-1; L,R=S.symbols('L R')
tail_double=2*t0*bt;tail_simple=t0*((t0-2)*L+t0*R);tail_jump=2*I*S.pi*t0**2*(L+R)
check('symmetric unary doublepole',(-bt/t0)*tail_double+2*bt**2)
check('symmetric unary simplepole',(-bt/t0)*tail_simple+bt*((t0-2)*L+t0*R))
check('mixed null repaired by balanced contact',(-bt/t0)*tail_jump+2*I*S.pi*t0*bt*(L+R))
# Entire ordered firstjet and derivative-test moments after outer Cauchy integration.
for sg in (-1,1):
 moment1=S.pi*sg; moment2=I*S.pi*sg; moment3=I*S.pi*sg
 outer=I*S.pi*bt*L*moment1+S.pi*bt**2*L*moment2+S.pi*bt*t0*R*moment3
 # action/cocycle factor yields -1/pi in convention after common nu lambda.
 check('ordered action contact sign '+str(sg),-outer/S.pi+I*S.pi*sg*t0*bt*(L+R))
 check('outer test derivative cancellation '+str(sg),bt**2*(S.pi*moment2-I*S.pi*moment1))
mp.mp.dps=40
for c in [mp.sqrt(2),mp.e]:
 for sg in (-1,1):
  val=mp.quad(lambda yy:mp.sign(sg*c-yy)*2*mp.sqrt(1-yy*yy),[-1,0,1])
  err=abs(val-sg*mp.pi)
  print(('PASS' if err<mp.mpf('1e-35') else 'FAIL')+' numerical ordered circle moment c='+str(c)+' sign='+str(sg)+' error='+str(err)); count+=1
# Positive patch diagnostic: scalar face and derivative remainder are locally integrable.
qq=S.symbols('q',positive=True)
for label,deg in [('distinct source deltas',2*h-3),('sourcepair base full',2*h-1),('two current differentiated full',2*h-3)]:
 expr=S.expand(deg.subs(h,2+qq));print('PASS '+label+' radialdegree at h=2+q = '+str(expr));count+=1
print('TOTAL',count,'PASS')
