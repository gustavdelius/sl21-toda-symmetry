"""Independent finite-part tests of triple angular and endpoint subtractions.

The radial test uses four independent nonconstant coefficient fields. It tests
the subtraction formula, not a full massive four-current correlation function.
"""
import mpmath as m
import sympy as S
m.mp.dps=45
passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    passed+=int(bool(ok)); failed+=int(not bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)

x=S.symbols('x')
rho=-m.sqrt(2)/10; t=2/rho; eta=t-2
P0=2*(x*x-x+1)
P1=(x-2)*(x+1)*(2*x-1)/3
P2=2*(x*x-x+1)**2/9
PS=4*(x*x-x+1)**2/3
polys={'I0':P0,'J0':P0*(1+x)**2/9,'I1r':P1*(1+x)/3,'IS2':PS,'IM2':P2}
def coeffs(poly):return [m.mpf(str(S.N(S.expand(poly).coeff(x,k),50))) for k in range(S.degree(poly,x)+1)]

# A direct angular integral after subtracting both endpoint jets, compared
# with beta-function continuation. Residual quadrature is genuinely convergent.
def half_fp(poly):
    cc=coeffs(poly)
    def h(v):return m.polyval(list(reversed(cc)),v)*(1-v)**(rho-2)
    h0=cc[0]; h1=(cc[1] if len(cc)>1 else 0)-(rho-2)*cc[0]
    series=[sum(cc[j]*(-1)**(k-j)*m.binomial(rho-2,k-j) for j in range(min(k,len(cc)-1)+1)) for k in range(18)]
    def rem(v):
        if v<m.mpf('.001'):
            return sum(series[k]*v**k for k in range(2,len(series)))
        return h(v)-h0-h1*v
    return m.quad(lambda v:v**(rho-2)*rem(v),[0,m.mpf('.001'),m.mpf('.5')])+h0*m.mpf('.5')**(rho-1)/(rho-1)+h1*m.mpf('.5')**rho/rho
def angular_beta(poly):
    return sum(c*m.beta(rho-1+k,rho-1) for k,c in enumerate(coeffs(poly)))
for name,poly in polys.items():
    actual=half_fp(poly)+half_fp(S.expand(poly.subs(x,1-x)))
    target=angular_beta(poly)
    check('continued angular finite part '+name,abs(actual-target)<m.mpf('1e-20'),
          'error='+m.nstr(abs(actual-target),4))

# The exact angular transform of polynomial times exp(v*x), obtained from
# Euler's integral in its convergence domain before continuation.
def transform(poly,v):
    return sum(c*m.beta(rho-1+k,rho-1)*m.hyp1f1(rho-1+k,2*rho-2+k,v)
               for k,c in enumerate(coeffs(poly)))
def integ_field(poly,k,r):
    # y in [r*(2-x)/3, 1-r*(1+x)/3], field exp(k*y).
    return (m.exp(k*(1-r/3))*transform(poly,-k*r/3)
            -m.exp(2*k*r/3)*transform(poly,-k*r/3))/k
k0,k1,ka,kb=map(m.mpf,('.23','-.17','.11','-.27'))
def A(r):
    return (-eta*(eta-1)*integ_field(P0,k0,r)-eta*r*integ_field(P1,k1,r)
            -eta*r*r*(integ_field(PS,ka,r)+integ_field(P2,kb,r)))
I0=angular_beta(P0); J0=angular_beta(polys['J0']); I1r=angular_beta(polys['I1r'])
IS=angular_beta(PS); IM=angular_beta(P2)
field_integral=lambda k:m.expm1(k)/k
A0=-eta*(eta-1)*I0*field_integral(k0)
A1=eta*(eta-1)*I0*(1+m.exp(k0))/2
A2=(-eta*(IS*field_integral(ka)+IM*field_integral(kb))
    -eta*(eta-1)*J0*k0*m.expm1(k0)/2+eta*I1r*m.expm1(k1))
for j,pred in enumerate((A0,A1,A2)):
    actual=m.diff(A,0,j)/m.factorial(j)
    check(f'triple radial coefficient from independent angular transform, level {j}',abs(actual-pred)<m.mpf('1e-32'))

# Evaluate the remainder via a convergent Taylor expansion, checked against
# the exact transform. This avoids catastrophic cancellation at tiny r.
series=m.taylor(A,0,10)
def rem(r):
    if r<m.mpf('.001'):return sum(series[j]*r**j for j in range(3,len(series)))
    return A(r)-A0-A1*r-A2*r*r
check('triple remainder evaluation agrees with exact transform',
      abs(sum(series[j]*m.mpf('.001')**j for j in range(3,len(series)))-(A(m.mpf('.001'))-A0-A1*m.mpf('.001')-A2*m.mpf('.001')**2))<m.mpf('1e-31'))
exponent=3*rho-3
target=m.quad(lambda r:r**exponent*rem(r),[0,m.mpf('.001'),1])+sum(aa/(3*rho-2+j) for j,aa in enumerate((A0,A1,A2)))
errors=[]
for eps in map(m.mpf,('.08','.04','.02','.01')):
    raw=m.quad(lambda r:r**exponent*A(r),[eps,1])
    ct=sum(aa*eps**(3*rho-2+j)/(3*rho-2+j) for j,aa in enumerate((A0,A1,A2)))
    err=abs(raw+ct-target); errors.append(err)
    print('DATA triple epsilon='+str(eps)+' finite-part='+m.nstr(raw+ct,14)+' error='+m.nstr(err,7),flush=True)
check('triple including descendants and endpoints converges',all(a>b for a,b in zip(errors,errors[1:])))
check('triple remainder has degree 3*rho+1',abs(errors[-1]/errors[-2]-2**(-(3*rho+1)))<m.mpf('.01'))
print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
