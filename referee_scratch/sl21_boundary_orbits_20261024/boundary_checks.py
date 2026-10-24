"""Endpoint covariance and anchored canonical weights from field data.

These are local or conditional results. There is no assumed charge algebra,
vacuum condensate, compactification lattice, or physical boundary limit.
"""
import sympy as S

count=0
def check(name,value):
    global count
    if not bool(value):
        print('FAIL',name)
        raise AssertionError(name)
    count+=1
    print('PASS',name)
def eq(x,y=0): return S.simplify(x-y)==0

b=S.symbols('b',positive=True)
t=b*b
a0=S.Matrix([S.I*b/S.sqrt(2),b/S.sqrt(2),S.I])
a2=S.Matrix([S.I*b/S.sqrt(2),-b/S.sqrt(2),-S.I])
T=a0+a2; d=-T/t; k=d+T; v=d+a2
R=S.Matrix([0,S.sqrt(2),S.I*b])
Sd=a0-a2
for name,p in [('Y0',a0),('Y2',a2),('V1',-T),('P',T),
               ('JA',-d),('JB exponential',d),('JB derivative',a0),
               ('K left',d+a2),('K right',a2),('L left',k),
               ('L right',T),('U',v),('F',k),('FA',-k)]:
    check('free-current neutrality '+name,eq(R.dot(p)))
check('orthogonal coordinates and nondegeneracy exceptions',
      eq(T.dot(Sd)) and eq(T.dot(R)) and eq(Sd.dot(R))
      and eq(T.dot(T),-2*t) and eq(Sd.dot(Sd),2*(t-2))
      and eq(R.dot(R),2-t))

n0,n2=S.symbols('n0 n2',integer=True)
N=n0+n2
shift=S.Matrix([S.sqrt(2)*S.pi*N/b,0,S.pi*(n0-n2)])
check('screening period lattice representative',
      eq(a0.dot(shift),2*S.pi*S.I*n0) and eq(a2.dot(shift),2*S.pi*S.I*n2))
check('even-root and composite period',eq(T.dot(shift),2*S.pi*S.I*N))
check('dual momentum period',eq(d.dot(shift),-2*S.pi*S.I*N/t))
# The differences are integer multiples of 2 pi i. This avoids treating a
# chosen logarithm of s as a globally single-valued physical operator.
for name,p,character,integer in [('U',v,1,n2),('F',k,1,N),('FA',-k,-1,-N)]:
    check(name+' period character s^N with appropriate sign',
          eq(p.dot(shift)+2*S.pi*S.I*character*N/t,2*S.pi*S.I*integer))
check('U and F have equal character but distinct UV momenta',
      v!=k and all(eq(x) for x in v-k+a0))
unit=shift.subs({n0:1,n2:0})
p0=a0.copy();p0[2]=0
p2=a2.copy();p2[2]=0
check('elementary period needs the chiral fermion flip',
      eq(p0.dot(unit),S.pi*S.I) and eq(p2.dot(unit),S.pi*S.I)
      and eq(unit[2],S.pi))

# Canonical derivation: integrate delta over (-infinity,x), for L<x<R.
# [Phi(y),phi(x)] = -2 pi i theta(x-y), opposite for bar phi.
for q in [-1,1]:
    left_phi=-2*S.pi*S.I
    right_phi=0
    wl=left_phi*T.dot(q*d)/(2*S.pi*S.I)
    wr=right_phi*T.dot(q*d)/(2*S.pi*S.I)
    check(f'q={q} anchored endpoint weights',eq(wl,-2*q) and eq(wr) and eq(wr-wl,2*q))
    check(f'q={q} antichiral weights reverse',eq(-wl,2*q))
    for name,pl,pr in [('J',q*d,S.zeros(3,1)),
                       ('Y tail',q*d+a2,a2),('two-source tail',q*d+T,T),
                       ('even tail',q*d-T,-T)]:
        check(f'q={q} {name} has unchanged disorder weight',
              all(eq(x) for x in pl-pr-q*d))
    for name,p,sgn,integer in [('U',v,1,1),('F',k,1,2),('FA',-k,-1,-2)]:
        # Endpoint-left BCH exponent is -2 pi i p.qd.
        exponent=-2*S.pi*S.I*p.dot(q*d)
        desired=4*S.pi*S.I*q*sgn/t # s^(-2q*sgn)
        check(f'q={q} {name} left endpoint phase',
              eq(exponent-desired,-2*S.pi*S.I*q*integer))
        check(f'q={q} {name} right endpoint commutes',eq(right_phi*p.dot(q*d)))

# Finite sector-shift model checks weight/phase signs without imposing any
# Serre or mixed algebra. Zero at truncation edges does not affect identities.
s=S.symbols('s',nonzero=True)
labels=list(range(-4,5))
left=S.diag(*labels)
right=S.zeros(len(labels))
top=right-left
EL=S.diag(*[s**n for n in labels])
ER=S.eye(len(labels))
for q in [-1,1]:
    Q=S.zeros(len(labels))
    for j,n in enumerate(labels):
        if n-2*q in labels: Q[labels.index(n-2*q),j]=1
    check(f'sector shift q={q} topological weight',top*Q-Q*top==2*q*Q)
    check(f'sector shift q={q} left endpoint phase',
          all(eq(x) for x in EL*Q-s**(-2*q)*Q*EL))
    check(f'sector shift q={q} right endpoint fixed',ER*Q==Q*ER)

nr,nt=S.symbols('nr nt',integer=True)
check('right-fixed BB boundary difference',
      eq(s**nr-s**(nr-nt),s**nr*(1-s**(-nt))))
check('right-fixed AA boundary difference',
      eq(s**(-nr)-s**(-nr+nt),s**(-nr)*(1-s**nt)))
nu=S.symbols('nu',nonzero=True)
for q in [-1,1]:
    check(f'K=nu^N weight q={q}',eq(nu**(2*q),(nu**2)**q))

dimF=-k.dot(k); dimU=-v.dot(v)
dimV=-T.dot(T)
dimlam1=2-dimV
check('F and U free-UV dimensions',eq(dimF,2*(t-1)**2/t) and eq(dimU,2/t-1))
check('AA endpoint dimension',eq(dimlam1+dimF,2/t-2))
check('BB endpoint dimension',eq(1+dimU,2/t))
check('scalar BB subtraction dimension',eq(2-(2*t-2)+dimF,2/t))
check('derivative BB subtraction dimension',eq(2-2*t+2+dimF,2/t))
check('endpoint dimensions agree with twice candidate charge spins',
      eq(dimlam1+dimF,2*(-d.dot(d)/2-1)) and eq(1+dimU,-d.dot(d)))

# Reconstruct the real-contour classical body potential from the roots.
u,mu=S.symbols('u mu',real=True,positive=True)
Vbody=mu**2*(S.exp(S.sqrt(2)*u)+S.exp(-S.sqrt(2)*u))
check('body potential derivative is a sinh with one real zero',
      eq(S.diff(Vbody,u),(2*S.sqrt(2)*mu**2*S.sinh(S.sqrt(2)*u)).rewrite(S.exp)))
check('body vacuum Hessian is positive',S.diff(Vbody,u,2).subs(u,0)==4*mu**2)
gc=S.symbols('gc',positive=True)
mapping={b:S.I*gc/S.sqrt(4*S.pi)}
Phi=S.Matrix([S.sqrt(4*S.pi)*u/gc,0,0])
check('source-to-canonical topological orientation',eq(T.subs(mapping).dot(Phi),-S.sqrt(2)*u))
check('complex kink 0 to i sqrt(2) pi has topological label -1',
      eq((-S.sqrt(2)*(S.I*S.sqrt(2)*S.pi))/(2*S.pi*S.I),-1))
print(f'{count} PASS, 0 FAIL')
