"""Independent reduction of Fateev 2001 (131)-(132), rank one."""
import sympy as s
import mpmath as mp

passed=0
def check(label,condition):
    global passed
    if not condition: raise AssertionError('FAIL '+label)
    passed+=1
    print('PASS '+label)
def zero(label,x): check(label,s.simplify(x)==0)

z,e=s.symbols('z e', positive=True)
v=z*z
alpha=s.sqrt(1+v)
Q=1/(2*z)
a=z+1/z+e
ap=2*Q-a
x=a-Q
xp=ap-Q
eta=-alpha
zero('reflected momentum',ap+z+e)
arguments=lambda xx:[1-2*xx*z,1-xx/z,s.Rational(1,2)-xx*z+alpha*eta,s.Rational(1,2)-xx*z-alpha*eta]
expected_a=[-2*v-2*z*e,-1/(2*v)-e/z,-1-2*v-z*e,1-z*e]
expected_ap=[2+2*v+2*z*e,2+1/(2*v)+e/z,z*e,2+2*v+z*e]
for i,(got,wanted) in enumerate(zip(arguments(x),expected_a)):
    zero('original gamma argument '+str(i),got-wanted)
for i,(got,wanted) in enumerate(zip(arguments(xp),expected_ap)):
    zero('reflected gamma argument '+str(i),got-wanted)
zero('source power uses barred momentum',xp/z-x/z+2+1/v+2*e/z)

# OPE scalar product uses independent physical kinetic normalization.
c=s.Matrix([-z-e,-s.I*alpha])
screen=s.Matrix([z,s.I*alpha])
zero('reflected/source OPE pairing',c.dot(screen)-(1-z*e))
zero('fused momentum radial component',(c+screen)[0]+e)
zero('fused momentum angular component',(c+screen)[1])
zero('radial integral pole denominator',2-2*c.dot(screen)-2*z*e)

# Generic-coupling resonance equations; descendant level is nonnegative.
n,m,D,ell=s.symbols('n m D ell', integer=True)
Ar=n-1-2*m
ref_const=D**2+n+2*m+ell-1
ref_v=D**2-Ar**2+4*m
zero('reflected dimension polynomial',
     (1+v)*D**2-v*Ar**2+n+(2+4*v)*m+ell-1-(ref_const+v*ref_v))
Au=1+n-2*m
u_const=D**2-1-n+6*m+ell
u_v=D**2-Au**2+4*m
zero('original dimension polynomial',
     (1+v)*D**2-v*Au**2-2*Au-1/v+n+(2+4*v)*m+ell-(-1-1/v)-(u_const+v*u_v))
zero('original resonance elimination',u_v-u_const-(Au*(1-Au)-ell))
solutions_ref=[]
solutions_u=[]
for n0 in range(9):
 for n2 in range(9):
  for mm in range(9):
   for level in range(5):
    subs={n:n0+n2,m:mm,D:n0-n2-1,ell:level}
    if ref_const.subs(subs)==0 and ref_v.subs(subs)==0: solutions_ref.append((n0,n2,mm,level))
    if u_const.subs(subs)==0 and u_v.subs(subs)==0: solutions_u.append((n0,n2,mm,level))
check('reflected resonance enumeration',solutions_ref==[(0,0,0,0),(1,0,0,0)])
check('original resonance enumeration',solutions_u==[(0,0,0,0)])

mp.mp.dps=80
for vv in [mp.sqrt(2),mp.e/2,mp.sqrt(3)/5]:
 zz=mp.sqrt(vv); aa=mp.sqrt(1+vv); qq=1/(2*zz); C=mp.sqrt(3)
 astar=zz+1/zz; et=-aa
 def amp(av):
  xx=av-qq
  return C**(xx/zz)*mp.gamma(1-2*xx*zz)*mp.gamma(1-xx/zz)/(mp.gamma(mp.mpf('.5')-xx*zz+aa*et)*mp.gamma(mp.mpf('.5')-xx*zz-aa*et))
 denom=C**(1+1/(2*vv))*(-1-2*vv)*mp.gamma(-1/(2*vv))
 check('generic original finite amplitude '+mp.nstr(vv,8),abs(amp(astar)/denom-1)<mp.mpf('1e-65'))
 slope=zz*C**(-2-1/vv)*mp.gamma(2+1/(2*vv))/((-1-2*vv)*mp.gamma(-1/(2*vv)))
 errors=[]
 for ee in [mp.mpf('1e-5'),mp.mpf('1e-7'),mp.mpf('1e-9')]:
  rr=amp(2*qq-(astar+ee))/amp(astar+ee)
  errors.append(abs(rr/(ee*slope)-1))
 check('simple-zero scaling '+mp.nstr(vv,8),errors[-1]<mp.mpf('2e-7') and errors[1]<errors[0]/90 and errors[2]<errors[1]/90)
 mu=2*vv*C/mp.pi
 expected=C**(-1-1/vv)*mp.gamma(1+1/(2*vv))/mp.gamma(-1/(2*vv))
 zero_pole=rr*(-mp.pi*mu/(zz*ee))
 check('zero times identity pole finite '+mp.nstr(vv,8),abs(zero_pole/expected-1)<mp.mpf('2e-7'))
 print('DATA v='+mp.nstr(vv,15)+' slope='+mp.nstr(slope,15)+' candidate finite='+mp.nstr(expected,15))
print('TOTAL '+str(passed)+' PASS, 0 FAIL')
