import sympy as S
n=0
def check(name,e):
 global n
 ok=S.simplify(e)==0
 print(('PASS ' if ok else 'FAIL ')+name)
 if ok:n+=1
 else:print(S.factor(e))
t,h,r=S.symbols('t h r',nonzero=True);rho=2/t;b=t-1
G=S.Matrix([[-rho,1],[1,-1]]);d=S.Matrix([1,0]);a0=S.Matrix([0,1]);a2=-t*d-a0
sd=lambda a,c:(a.T*G*c)[0]
check('distinct unary fusion holomorphic exponent',-sd(d+a0,a2)+(h-t)-(h-2))
check('distinct unary fusion antiholomorphic exponent',-sd(a0,d+a2)+(h-t)-(h-2))
check('distinct-source two-face remaining denominator sign',1/(r*(-r))+1/r**2)
check('equal source unary fusion oscillator is regular',sd(d+a2,a2))
check('equal source factorial cancels two assignments',S.Rational(1,2)*(b*b+b*b)-b*b)
R,K,Rp,Kp=S.symbols('R K Rp Kp')
actual_simple=-R*Kp+b*Rp*K
check('spectator Laurent regulator coefficient is beta+1',actual_simple-(-(Rp*K+R*Kp)+t*Rp*K))
print(f'TOTAL {n} PASS')
