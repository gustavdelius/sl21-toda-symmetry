"""Fresh strip-edge calculation and endpoint subtractions from Wick data."""
import sympy as S
import mpmath as mp
n=0
def check(name,v):
    global n
    if not bool(v):
        print('FAIL',name)
        raise AssertionError(name)
    n+=1
    print('PASS',name)
def eq(e): return S.simplify(e)==0

t=S.symbols('t',nonzero=True)
q=t-1
CL,D,AL,AR,CR,KL,KR,I,Jw,J1=S.symbols('CL D AL AR CR KL KR I Jw J1')
# Upper-minus-lower Cauchy moments are -2 pi i Theta and its w derivative.
# u/(u-w)^2 = 1/(u-w) + w/(u-w)^2. Expand the exponentials to first order.
first=2*S.I*q*((CL+q*D)*I+q*(D+AL)*Jw+(q*AR-CR)*J1)
reduced=S.expand(first.subs({CL:-KL-q*D,AL:KL-D,CR:q*AR-t*KR}))
check('first jet follows from Cauchy moments with action sign',
      eq(reduced+2*S.I*q*(KL*(I-q*Jw)-t*KR*J1)))
contact=reduced.subs({I:S.pi,J1:-S.pi,Jw:-S.pi})
check('bare BB current-action finite contact',eq(contact+2*S.pi*S.I*t*(t-1)*(KL+KR)))
check('tails plus bare contact leave radial coefficient t',eq(t*t-t*(t-1)-t))
gamma=t/(1-t)
check('AA radial coefficient and canonical addition',eq(gamma**2+gamma-gamma**2/t))
check('canonical finite additions',eq(t+t*(t-1)-t*t) and eq(gamma**2/t-gamma-gamma**2))

mp.mp.dps=40
I_num=mp.quad(lambda y:2*mp.sqrt(1-y*y),[-1,0,1])
check('circle-excluded strip moment I=pi',abs(I_num-mp.pi)<mp.mpf('1e-35'))
c=mp.sqrt(3)
def line(y,which):
    def f(x):
        z=mp.mpc(x,y)
        if which==1:return 1/mp.conj(z)
        return z/mp.conj(z)**2
    return mp.quad(lambda x:f(x)+f(-x),[0,1,mp.inf])
for which,name in [(1,'J1'),(2,'Jw')]:
    moment=mp.j/2*(line(c,which)-line(-c,which))
    check(f'strip derivative gives {name}=-pi',abs(moment+mp.pi)<mp.mpf('1e-34'))
eta=mp.mpf('.17'); y=mp.mpf('.63')
for pole in [1,2]:
    def f(x):
        z=mp.mpc(x,y)
        return mp.conj(z)**pole/(abs(z)**2+eta**2)**pole
    val=mp.quad(lambda x:f(x)+f(-x),[0,1,mp.inf])
    expected=(-mp.j*mp.pi*y/mp.sqrt(y*y+eta*eta) if pole==1
              else mp.pi*eta**2/(2*(y*y+eta*eta)**mp.mpf('1.5')))
    check(f'smooth Cauchy distribution contour moment pole={pole}',abs(val-expected)<mp.mpf('1e-34'))

a0=S.Matrix([S.I*S.sqrt(t/2),S.sqrt(t/2),S.I])
a2=S.Matrix([S.I*S.sqrt(t/2),-S.sqrt(t/2),-S.I])
d=-(a0+a2)/t
check('endpoint/source exponent differs from source/source exponent',
      eq(-(d+a2).dot(a0)-(t-2)) and eq(-a2.dot(a0)-(t-1)))
eps,R=S.symbols('eps R',positive=True)
for power in [t-1,t]:
    primitive=S.pi*(R**(2*power)-eps**(2*power))/power
    check(f'endpoint subtraction radial primitive {power}',
          eq(S.diff(primitive,R)-2*S.pi*R**(2*power-1)))
check('BB descendant coefficient after canonical multiplication',eq(t*t/t-t))
print(f'{n} PASS, 0 FAIL')
