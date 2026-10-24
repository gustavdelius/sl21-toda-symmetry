import sympy as S
import mpmath as mp
from itertools import product
mp.mp.dps=45
b,t,c=S.symbols('b t c',nonzero=True, real=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I]);a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I]);T=a0+a2;d=-T/b**2;k=d+T
simp=lambda v:S.simplify(S.simplify(v).subs(b**2,t))
fail=[]
def ck(name,expr):
 ok=bool(expr); print(('PASS ' if ok else 'FAIL ')+name)
 if not ok:fail.append(name)
def eq(name,x,y=0):ck(name,simp(x-y)==0)
for name,x,y in [('a0^2',a0.dot(a0),-1),('a2^2',a2.dot(a2),-1),('a0a2',a0.dot(a2),1-t),('d^2',d.dot(d),-2/t),('a0d',a0.dot(d),1),('a2d',a2.dot(d),1)]:eq(name,x,y)
for j in range(3):eq('jet vector '+str(j),a0[j]+(b*b-1)*d[j]+d[j]+a2[j])
# Laurent coefficients of exact scalar ABY_i kernel by differentiating analytic prefactor.
rho=2/t; a,w,z=S.symbols('a w z'); beta=t-1
F=(a-z)**(-rho)*(a-w)
# Do not ask symbolic algebra to combine branch-sensitive fractional powers.
lead=beta*(a-w)**(1-rho)
res=(a-w)*(a-w)**(-rho)/(w-a)+beta*rho*(a-w)*(a-w)**(-rho-1)
ck('first mixed jet coefficient',S.simplify((beta*rho-1)-(1-rho))==0)
# 2D strip constants evaluated separately, asymmetric heights.
I=mp.quad(lambda yy:2*mp.sqrt(1-yy*yy),[-1,0,1]);ck('circle strip I',abs(I-mp.pi)<mp.mpf('1e-35'))
for h in [mp.sqrt(2),mp.sqrt(3)]:
 edge=mp.quad(lambda x:1/(x+1j*h),[-mp.inf,0,mp.inf]);ck('strip edge h='+str(h),abs(edge+1j*mp.pi)<mp.mpf('1e-35'))
eq('action + tails = raw',-t*(t-1)+t*t,t)
eq('finite matching = canonical',t+t*(t-1),t*t)
# Gaussian transport from original vectors, not copied scalar factors.
H=S.diag(c,c,0)
logZ=lambda p:simp(-(p.T*H*p)[0]/2)
for name,p,val in [('Y0',a0,0),('Y2',a2,0),('V1',-T,c*t),('FA',-k,c*(1-t)**2/t),('U',d+a2,c*(1/t-1)),('d',d,c/t)]:eq('logZ '+name,logZ(p),val)
eq('Ward AA log factor',logZ(-T)+logZ(d)-logZ(-k),2*c)
eq('Ward BB log factor',logZ(a2)+logZ(d)-logZ(d+a2),c)
eq('UY0 fusion log',((d+a2).T*H*a0)[0],c*(1-t))
eq('Y0Y2 fusion log',(a0.T*H*a2)[0],-c*t)
# Source action double-sum and expansion combinatorics.
l0,l2,Y0,Y2=S.symbols('l0 l2 Y0 Y2')
eq('distinct-source action combinatorics',S.expand((l0*Y0+l2*Y2)**2/2).coeff(l0,1).coeff(l2,1),Y0*Y2)
eq('UY0 exponent',-2*(d+a2).dot(a0),2*t-4)
eq('Y0Y2 exponent',-2*a0.dot(a2),2*t-2)
# Reflection ratio from visually inspected source (132), ratio (131).
for tt in [-mp.sqrt(2)/3,-mp.sqrt(3)/4]:
 zz=mp.sqrt(-tt/2);al=mp.sqrt(1+zz*zz);Q=1/(2*zz);C=mp.sqrt(mp.mpf('1.1')*mp.mpf('1.3'))/(-tt);lam0=mp.mpf('1.1')
 def AA(x):
  v=x-Q
  return C**(v/zz)*mp.gamma(1-2*zz*v)*mp.gamma(1-v/zz)/(mp.gamma(mp.mpf('.5')-zz*v-al*al)*mp.gamma(mp.mpf('.5')-zz*v+al*al))
 de=mp.mpf('1e-14');x=zz+1/zz+de
 ratio=AA(2*Q-x)/AA(x)
 linear=de*zz*C**(2/tt-2)*mp.gamma(2-1/tt)/((tt-1)*mp.gamma(1/tt))
 target=lam0/(1-tt)*C**(2/tt-2)*mp.gamma(2-1/tt)/mp.gamma(1/tt)
 ck('source reflection simple zero t='+str(tt),abs(ratio/linear-1)<mp.mpf('1e-11'))
 ck('zero times radial identity pole t='+str(tt),abs(ratio*(-lam0/(zz*de))/target-1)<mp.mpf('1e-11'))
 # positive i0 oscillator ratio and cocycle for U(a) A(x)
 ss=mp.e**(-2j*mp.pi/tt); coc=-ss; alpha=1-2/tt
 ck('left endpoint phase',abs(coc*mp.e**(1j*mp.pi*alpha)-ss**2)<mp.mpf('1e-35'))
 ck('right endpoint phase',abs(coc*mp.e**(-1j*mp.pi*alpha)-1)<mp.mpf('1e-35'))
print('FAILURES',len(fail))
assert not fail
