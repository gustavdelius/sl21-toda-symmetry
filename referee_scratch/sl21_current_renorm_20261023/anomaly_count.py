"""Kinematic anomaly count in mass-analytic Gaussian conformal perturbation theory.

This is not a proof of a global charge or of convergence of the expansion.
The algebraic argument in the report accompanies the finite enumeration.
"""
import sympy as s


def check(name,value):
    ok=value if isinstance(value,bool) else s.simplify(value)==0
    print(('PASS' if ok else 'FAIL')+' '+name,flush=True)
    assert ok,name


a,c,ne,t=s.symbols('a c ne t')
n0,n2=a+ne,c+ne
mass=n0+n2+2*ne
ferm=n0-n2
S=a+c
basebar=t*a*c+ferm**2/2
basehol=1/t-S+basebar
coupling=(n0+n2)/2+ne*(1-t)
lh=s.expand(1+1/t-coupling-basehol)
lb=s.expand(1-coupling-basebar)
check('left-right descendant difference',lh-lb-S)
check('common coefficient of t',s.diff(lh,t)-(ne-a*c))
check('antiholomorphic generic-coupling integer condition',
      lb.subs(ne,a*c)-(1-S*(S+1)/2))
check('holomorphic generic-coupling integer condition',
      lh.subs(ne,a*c)-(1-S*(S-1)/2))

allowedS=[n for n in range(-20,21)
          if 1-n*(n+1)//2>=0 and 1-n*(n-1)//2>=0]
check('nonnegative levels require S=-1,0,1',allowedS==[-1,0,1])
solutions=[]
for aa in range(-20,21):
    for cc in range(-20,21):
        e=aa*cc
        m0=aa+e
        m2=cc+e
        if e>=0 and m0>=0 and m2>=0 and aa+cc in allowedS:
            solutions.append((m0,m2,e))
check('generic integer solutions contain only vacuum and single Yukawas',
      sorted(solutions)==[(0,0,0),(0,1,0),(1,0,0)])

print('ASSUMPTIONS: analytic coupling expansion; UV Fock fields; no independent finite P coupling',flush=True)
print('RESULT: ordinary local higher-order anomalies are kinematically excluded for generic irrational t',flush=True)
for m0,m2,e in [(3,0,0),(2,1,0),(1,2,0),(0,3,0),(1,0,1),(0,1,1)]:
    repl={a:m0-e,c:m2-e,ne:e}
    Lh=s.simplify(lh.subs(repl))
    Lb=s.simplify(lb.subs(repl))
    print(f'ORDER3 {(m0,m2,e)}: ({Lh}, {Lb})',flush=True)
check('first triple-Yukawa resonance at t=-1/2',
      lb.subs({a:2,c:1,ne:0,t:-s.Rational(1,2)}))
check('first mixed even/Yukawa resonance at t=2 has level (0,1)',
      lh.subs({a:0,c:-1,ne:1,t:2})**2
      +(lb.subs({a:0,c:-1,ne:1,t:2})-1)**2)
