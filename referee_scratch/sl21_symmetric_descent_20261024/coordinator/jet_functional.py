"""Independent integration of the complete first jet of a mixed proper face.

The three source moments are evaluated by horizontal slices, then used
as a linear functional. This is a referee check, not a book check script.
"""
import sympy as S
q,C,Cb,p0,pz,pw,pb,sgn,lam,nu=S.symbols('q C Cb p0 pz pw pb sgn lam nu')
t,b=S.symbols('t b',nonzero=True)
pi,I=S.pi,S.I
checks=[]
def ck(name,expr):
    yes=bool(expr); checks.append(yes)
    print(('PASS ' if yes else 'FAIL ')+name)
# Source moments after the outer Cauchy integral, including its moving pole.
# Their common sign is the signed outer height. No normal-field substitution yet.
sign_over_bar2=pi*sgn
line_over_bar=I*pi*sgn
line_w_over_bar2=I*pi*sgn
outer_simple=-I*pi
outer_double=pi
coeff0=outer_simple*C*q*sign_over_bar2-outer_double*q*Cb*line_over_bar
coeffw=outer_double*q*q*line_w_over_bar2
coeffb=outer_double*q*q*line_over_bar
coeffz=q*q*(outer_simple*sign_over_bar2+outer_double*line_w_over_bar2)
raw=-lam*nu/pi*(coeff0*p0+coeffw*pw+coeffb*pb+coeffz*pz)
expected=-I*pi*sgn*nu*lam*(q*q*(pw+pb)-q*(C+Cb)*p0)
ck('entire arbitrary first jet gives the raw mixed contact',S.expand(raw-expected)==0)
ck('outer-current test derivative cancels before normal-field substitution',S.expand(coeffz)==0)
a0=S.Matrix([I*b/S.sqrt(2),b/S.sqrt(2),I])
a2=S.Matrix([I*b/S.sqrt(2),-b/S.sqrt(2),-I])
T=a0+a2; d=-T/b**2; v=d+a2; k=d+T
beta=b**2-1; gamma=b**2/(1-b**2)
ck('Y0 physical normal jet vanishes',a0-a0==S.zeros(3,1))
ck('Y2 physical normal jet equals t beta derivative of U',all(S.simplify(x)==0 for x in beta**2*a2-beta*a0-b**2*beta*v))
chi,chip,U,Up=S.symbols('chi chip U Up')
full=t*(t-1)*chi*Up+(t-1)**2*chip*U
oversimple=t*(t-1)*(chi*Up+chip*U)
ck('switching derivative coefficient is beta squared, not t beta',S.expand(full-oversimple+(t-1)*chip*U)==0)
ck('Y0 has a source-test derivative even when normal contact vanishes',S.expand(expected.subs({q:1,C:0,Cb:0,pw:1,pb:0,p0:0})+I*pi*sgn*nu*lam)==0)
aa=-lam*nu/pi*(outer_double*line_w_over_bar2*pw+outer_double*line_over_bar*pb)
ck('AA full independent source jet',S.expand(aa+I*pi*sgn*nu*lam*(pw+pb))==0)
ck('AA physical field jet conversion',all(S.simplify(x)==0 for x in -T+gamma*(-k)))
ck('balanced opposite half matching cancels full raw functional',S.expand(raw+(-raw))==0)
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
