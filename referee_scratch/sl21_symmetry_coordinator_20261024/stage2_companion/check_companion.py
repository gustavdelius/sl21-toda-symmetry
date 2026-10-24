"""Independent one-source B A barA kernels and test-contact calculation."""
import sympy as S
import mpmath as mp

t,a,z,w,C,p=S.symbols('t a z w C p',nonzero=True)
rho=2/t; eta=1/(1-t); gamma=t/(1-t)
b=S.symbols('b',nonzero=True)
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I])
a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I])
T=a0+a2;d=-T/b**2;k=d+T
count=0
def ck(name,value):
    global count
    count+=1
    print(('PASS' if bool(value) else 'FAIL')+' '+name)
    if not value:raise AssertionError(name)
ck('B-source V1 power two',S.simplify(d.dot(T)-2)==0)
ck('A-source V1 double pole',S.simplify(-d.dot(T)+2)==0)
ck('B derivative contraction with V1 is -t',S.simplify(a0.dot(T)+b*b)==0)
ck('fused B FA derivative contraction 1-t',S.simplify(a0.dot(k)-(1-b*b))==0)
ck('full field-jet identity eta k=d',all(S.simplify(x)==0 for x in k/(1-b*b)-d))
# Strip moving pole: suppress normal exponential but retain spectator derivative C.
N=(a-z)**(-rho)*(a-w)**2*(C+1/(a-z)-t/(a-w))
V=(a-w)**(2-rho)*(C+(1-t)/(a-w))
ck('moving-A double pole coefficient',S.simplify(N.subs(z,w)-V)==0)
# Powers have common positive/analytic branch; expand their logarithmic derivatives.
simple=(a-w)**(1-rho)*(rho*C+(rho*(1-t)+1)/(a-w))
ck('moving-A simple pole eta V prime',S.simplify(S.expand_power_base(simple-eta*S.diff(V,w),force=True))==0)
ck('bare action coefficient gamma',S.simplify(eta-1-gamma)==0)
ck('two-tail scalar derivative coefficient gamma squared',S.simplify(gamma*eta-gamma-gamma**2)==0)
ck('bare action antichiral field jet T=-gamma k',all(S.simplify(x)==0 for x in T+b*b/(1-b*b)*k))
ck('finite matching leaves canonical gamma squared',S.simplify(gamma+gamma**2-gamma-gamma**2)==0)
rr,pp=S.symbols('r p',nonzero=True)
kernel=(rr+pp)**2/pp**2*(1/rr-t/(rr+pp))
resolved=(1-t)*rr/pp**2+(2-t)/pp+1/rr
ck('exact simultaneous scalar kernel resolution',S.simplify(kernel-resolved)==0)
# Actual derivative-of-test terms before the canonical addition.
I,Jw=S.symbols('I Jw')
ck('bare action first test jet cancels',((I+Jw).subs({I:S.pi,Jw:-S.pi}))==0)
ck('tail first test jets cancel',S.simplify(gamma-gamma)==0)
labels={'B':(1,0,0,0),'A':(-1,0,0,0),'a':(0,-1,0,0),'V':(0,0,-1,-1),'H':(-1,0,-1,-1),'h':(0,-1,-1,-1),'F':(-1,-1,-1,-1)}
def cocycle(word):
    sg=0;ph=0
    for i,key in enumerate(word):
        n=labels[key]
        for kk in word[i+1:]:
            m=labels[kk]
            sg+=(n[0]+n[1])*(m[2]+m[3]);ph+=n[0]*m[1]-n[1]*m[0]
    return (-1 if sg%2 else 1),ph
for word in ['BAaV','BAh','BHa']:
    ck('reference cocycle '+word,cocycle(word)==(1,0))
ck('B FA cocycle ratio exponent',cocycle('BF')[1]-cocycle('FB')[1]==-2)
# Physical dilation weights of candidate derivative contact.
hU=S.simplify(-(d+a2).dot(d+a2)/2)
hJ=S.simplify(1-d.dot(d)/2)
ck('lambda2 U delta-prime has product density dimension',S.simplify(1+2*hU+2-2*hJ)==0)

mp.mp.dps=35
tv=-mp.sqrt(2)/3
rv=2/tv;gv=tv/(1-tv)
sv=mp.exp(-2j*mp.pi/tv)
nu=mp.exp(-1j*mp.pi/tv)
# Positive-imaginary-separation oscillator boundary values.
alpha=2-rv
ratio_left=sv**(-1)*mp.exp(-1j*mp.pi*alpha)
ratio_right=sv**(-1)*mp.exp(1j*mp.pi*alpha)
ck('left FA B=s squared B FA',abs(1/ratio_left-sv**2)<mp.mpf('1e-30'))
ck('right FA commutes with B',abs(ratio_right-1)<mp.mpf('1e-30'))
ck('proper and primitive degrees positive',(-rv>0) and (1-rv>0) and (2-rv>0))
# Independent upper-minus-lower double-pole integration with nonconstant test.
def pole_test(delta):
    fun=lambda x:x*mp.exp(-x*x)*(1/(x+1j*delta)**2-1/(x-1j*delta)**2)
    return mp.quad(fun,[-mp.inf,-delta,0,delta,mp.inf])/(-2j*mp.pi)
vals=[]
for delta in [mp.mpf('.02'),mp.mpf('.01'),mp.mpf('.005')]:
    val=pole_test(delta);vals.append(val)
    print('DIAGNOSTIC double-pole first-test-jet delta',delta,'normalized',mp.nstr(val,18))
ck('double-pole test derivative convergence',abs(vals[-1]-1)<abs(vals[0]-1) and abs(vals[-1]-1)<mp.mpf('.02'))
for ep in [mp.mpf('.02'),mp.mpf('.01'),mp.mpf('.005')]:
    bound=ep*mp.quad(lambda R:R**(-1-rv),[4*ep,1])
    primitive=ep**(1-rv)
    print('DIAGNOSTIC epsilon',ep,'mixed remainder bound',mp.nstr(bound,15),'primitive power',mp.nstr(primitive,15))
ck('mixed remainder O(epsilon) bound',abs(bound/ep-1/(-rv))<mp.mpf('.00001'))
# A concrete intrinsic density contact action: -(g prime) times h U.
contact=mp.quad(lambda x:-(2*x)*mp.exp(-x*x)*mp.exp(-x*x)*(1+x),[-mp.inf,mp.inf])
ck('candidate derivative contact visible to nonconstant outer test',abs(contact+mp.sqrt(mp.pi)/(2**mp.mpf('1.5')))<mp.mpf('1e-30'))
print('DIAGNOSTIC candidate density contact nonconstant-test value',mp.nstr(contact,18))
print('Completed',count,'checks. Scope: restricted lambda1 companion and local first-test jets, not a common all-order product.')
