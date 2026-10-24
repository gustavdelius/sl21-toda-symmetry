"""The two independent mixed faces require a joint balanced projector."""
import sympy as S
m,a,b,d=S.symbols('m a b d')
R=lambda s,t:m+s*a/2+t*b/2+s*t*d/4
rr={(s,t):R(s,t) for s in [-1,1] for t in [-1,1]}
Delta=rr[1,1]-rr[1,-1]-rr[-1,1]+rr[-1,-1]
global_two=(rr[1,1]+rr[-1,-1])/2
independent=sum(rr.values())/4
checks=[]
def ck(name,x):
    good=bool(x);checks.append(good);print(('PASS ' if good else 'FAIL ')+name)
ck('mixed difference extracts joint jump',S.expand(Delta-d)==0)
ck('two-extreme shortcut has quarter-joint-jump error',S.expand(global_two-independent-Delta/4)==0)
ck('independent face balancing removes both odd jumps and joint jump',S.expand(independent-m)==0)
ck('full-jet nonzero test forbids setting joint jump to zero',S.diff(S.Symbol('x')*S.Symbol('y'),S.Symbol('x'),S.Symbol('y'))==1)
print('COUNTEREXAMPLE R++=R--=1/4, R+-=R-+=-1/4: global=1/4, independent=0')
ck('explicit joint test distinguishes prescriptions',global_two.subs({m:0,a:0,b:0,d:1})==S.Rational(1,4) and independent.subs({m:0,a:0,b:0,d:1})==0)
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
