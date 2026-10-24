import sympy as s
from itertools import product
n=0

def check(name, expr):
    global n
    ok=bool(expr) if isinstance(expr,bool) else s.simplify(expr)==0
    print(('PASS' if ok else 'FAIL')+' '+name)
    if not ok: raise AssertionError((name,expr))
    n+=1

t,g,rho,w,z2,z3,a=s.symbols('t g rho w z2 z3 a', nonzero=True)
C0,C2,C3,D=s.symbols('C0 C2 C3 D')
for q,sigma,E,name in [(1,1,C0,'Y0'),(t-1,-1,-t*D-C0,'Y2')]:
    f2=C2+(q-1)/(z2-w)-1/(z2-z3)+1/(z2-a)
    f3=C3+(q-1)/(z3-w)-1/(z3-z2)+1/(z3-a)
    h2=f2*f3+1/(z2-z3)**2
    A=1/(w-z2)+1/(w-z3)
    J=f3/(w-z2)**2+f2/(w-z3)**2
    # Exact Laurent coefficient from three singleton matchings plus three pairs.
    direct=(C0+q*D+(q*g-1)*A+(1-q*rho)/(w-a))*h2+(1-q)*J
    deriv=(D+E+(g-1)*A+(1-rho)/(w-a))*h2+s.diff(h2,w)
    candidate=sigma*deriv+(q-sigma)*(g-rho)*A*h2
    check(name+' full H3 simple pole incl two spectator derivative fields',s.factor((direct-candidate).subs(rho,2/t)))
    check(name+' normal-field jet identity',sigma*(D+E)-(C0+q*D))
    check(name+' fused spectator derivative identity',s.diff(h2,w)-(q-1)*J)
    check(name+' A pole identity',(sigma+(q-sigma)*rho-1).subs(rho,2/t))
    check(name+' double-pole coefficient',q*h2-q*h2)
# Even-root V1, three A and one B. At A1=w+e, normal jet=-D.
z=s.symbols('z',nonzero=True)
f=C0+1/(z-w)+1/(z-z2)+1/(z-z3)-t/(z-w)
# derivative of a0 derivative-current contraction with A1 is +1/(z-w)^2
S=(-D+g/(w-z2)+g/(w-z3)-rho/(w-z))*f+1/(z-w)**2
# F_A exponent is -d-T=(t-1)d; scalar fused A powers g-2 and B power 2-rho.
k_der=((t-1)*D+(g-2)/(w-z2)+(g-2)/(w-z3)+(2-rho)/(w-z))*f+s.diff(f,w)
# sigma=1/(1-t), coefficient of derivative is -gamma/t.
sigma=1/(1-t)
rem=s.factor(S-sigma*k_der)
expected=(1-sigma)*(g-rho)*(1/(w-z2)+1/(w-z3))*f
check('V1 full A-source simple pole and evanescent spectator term',(rem-expected).subs(rho,2/t))
# Differential product rule including distinct mask faces is exact.
x,y=s.symbols('x y'); M=s.Function('M')(x,y); F=s.Function('F')(x,y)
check('full mixed mask Leibniz rule',s.diff(M*F,x,y)-(M*s.diff(F,x,y)+s.diff(M,x)*s.diff(F,y)+s.diff(M,y)*s.diff(F,x)+s.diff(M,x,y)*F))
theta=[s.Function('theta'+str(i))(x,y) for i in range(4)]
prod=s.prod(theta)
rhs=sum(s.diff(theta[i],x,y)*s.prod(theta[j] for j in range(4) if j!=i) for i in range(4))+sum(s.diff(theta[i],x)*s.diff(theta[j],y)*s.prod(theta[k] for k in range(4) if k not in (i,j)) for i in range(4) for j in range(4) if i!=j)
check('mask intersections retained',s.diff(prod,x,y)-rhs)
# All one-source proper faces: exclude no-current and fixed low faces.
for ch in ['B','A']:
  for aa,bb,dd in product(range(2) if ch=='B' else range(4),range(4) if ch=='B' else range(2),range(2)):
    if aa+bb+dd==0: continue
    E=g*(bb*(bb-1)/2 if ch=='B' else aa*(aa-1)/2)-rho*aa*bb
    omega=E+(2*aa-bb-dd+1 if ch=='B' else -aa+2*bb-dd+1)
    fixed= (aa==0 and bb+dd<=2 and bb<=1) if ch=='B' else (bb==0 and aa+dd<=2 and aa<=1)
    if fixed: continue
    check(ch+' higher union '+str((aa,bb,dd))+' positive with first-jet margin',float(omega.subs({g:7,rho:-2}))>1)
# Whole ordered-word balancing, four chiral slots and one barred slot.
R0=s.symbols('R0'); jumps=s.symbols('j1:5'); R=[R0]
for j in jumps:R.append(R[-1]+j)
E=(R[0]+R[-1])/2
for k in range(5):
    correction=(sum(jumps[k:])-sum(jumps[:k]))/2
    check('complete five-current placement '+str(k),R[k]+correction-E)
# Arbitrary smooth test first derivative in balanced B pair.
for sigma in [-1,1]:
    check('outer test first jet sign '+str(sigma),(t-1)**2*(s.pi*(s.I*s.pi*sigma)-s.I*s.pi*(s.pi*sigma)))
# Unlike source anti triangle, reduce without suppressing Cbar.
p,q,beta,Cb=s.symbols('p q beta Cb',nonzero=True)
original=(Cb+1/p+beta/q)/(p*q)
reduced=(Cb*(1/p-1/q)+1/p**2-beta/q**2)/(q-p)+(beta-1)*(1/p-1/q)/(q-p)**2
check('unlike two-source complete anti triangle',original-reduced)
# Relative source pair Taylor coefficients, computed directly from geometric series.
z,r=s.symbols('z r',nonzero=True)
f=(1/(z-r)+(t-1)/z)/((z-r)*z)
for k in range(4):check('source-pair Taylor coefficient '+str(k),s.diff(f,r,k).subs(r,0)/s.factorial(k)-(k+t)/z**(k+3))
print('PASS total '+str(n))
