"""Fresh free-word checks. No imports of preceding investigation scripts.

The identities tested are necessary algebraic compatibility conditions, not
evaluations of renormalized three-contour integrals or a vacuum construction.
"""
from collections import defaultdict
import sympy as sp

s, v, a, b = sp.symbols('s v a b', nonzero=True)
c = 1+s+1/s
passed = 0

def check(label, value):
    global passed
    if isinstance(value, sp.MatrixBase):
        ok = all(sp.cancel(x) == 0 for x in value)
    elif isinstance(value, bool):
        ok = value
    else:
        ok = sp.cancel(value) == 0
    print(('PASS' if ok else 'FAIL')+' '+label)
    if not ok:
        raise AssertionError((label, value))
    passed += 1

# Endpoint brackets start from the anchored integral, with q_A=-1, q_B=+1.
# [n_L,Q_q]=-q T.d; [n_R,Q_q]=0. Recompute T.d from vectors.
bb = sp.symbols('bb', nonzero=True)
a0=sp.Matrix([sp.I*bb/sp.sqrt(2),bb/sp.sqrt(2),sp.I])
a2=sp.Matrix([sp.I*bb/sp.sqrt(2),-bb/sp.sqrt(2),-sp.I])
T=a0+a2
d=-T/bb**2
check('canonical T.d', T.dot(d)-2)
for name,q in [('A',-1),('B',1)]:
    check(name+' right-minus-left weight', q*T.dot(d)-2*q)

# Word defect: W Z - r_W Z W. For XZ=r_X ZX+f_X(K),
# each matching letter can contract. Move each endpoint K^m to the right.
weight={'A':-1,'B':1}
def cross_defect(poly, barred_label, endpoint_power):
    answer=defaultdict(lambda:sp.S.Zero)
    for word,coeff in poly.items():
        for i,letter in enumerate(word):
            if letter != barred_label:
                continue
            suffix=word[i+1:]
            braid_power=sum(1 if x==barred_label else -1 for x in suffix)
            endpoint_power_shift=endpoint_power*sum(weight[x] for x in suffix)
            remainder=word[:i]+word[i+1:]
            answer[remainder] += coeff*s**(braid_power+endpoint_power_shift)
    return {w:sp.factor(q) for w,q in answer.items()}

C={('B','A'):1,('A','B'):-1/s}
# f_A=a0+a2 K^2, f_B=b0+bminus2 K^-2; test each independent mode.
expected={('A',0):{'B':1-s**-2}, ('A',2):{'B':0},
          ('B',0):{'A':0}, ('B',-2):{'A':s-s**-1}}
for (bar,m),target in expected.items():
    actual=cross_defect(C,bar,m)
    for word,coef in actual.items():
        check('triple '+bar+' K^'+str(m)+' '+' '.join(word),coef-target[''.join(word)])

def serre(majority,minority):
    return {(majority,)*(3-j)+(minority,)+(majority,)*j:coef
            for j,coef in enumerate([1,-c,c,-1])}

polys={'B':serre('B','A'),'A':serre('A','B')}
for majority,poly in polys.items():
    for barred, modes in [('A',[0,2]),('B',[0,-2])]:
        for m in modes:
            for word,coef in cross_defect(poly,barred,m).items():
                check('quartic '+majority+' vs bar'+barred+' K^'+str(m)+' '+''.join(word),coef)

# Derive mode classification without installing the expected endpoint form.
y=sp.symbols('y',nonzero=True)
p=lambda x:1-c*x+c*x**2-x**3
check('single-minority polynomial factor',p(y)-(1-y)*(1-s*y)*(1-y/s))
check('majority coefficient factor',y+1+1/y-c-(y-s)*(y-1/s)/y)
for majority,poly in polys.items():
    for m in range(-4,5):
        out=cross_defect(poly,majority,m)
        vanishes=all(sp.cancel(x)==0 for x in out.values())
        permitted=(m in ([0,-2] if majority=='B' else [0,2]))
        check('majority mode classification '+majority+' m='+str(m), vanishes==permitted)

# Explicit associative 4x4 matrix witness. v is a symbolic indeterminate,
# avoiding accidental root-of-unity values; s=v^2. Matrix construction is
# independent of word reduction above and keeps both endpoint amplitudes free.
n=3
qnum=lambda j:(v**j-v**(-j))/(v-v**-1)
K=sp.diag(*[v**(2*j-n) for j in range(n+1)])
E=sp.zeros(n+1)
F=sp.zeros(n+1)
for j in range(n):
    E[j+1,j]=1
for j in range(1,n+1):
    F[j-1,j]=qnum(j)*qnum(n-j+1)
A=F
B=E
barA=a*(v-v**-1)/v**2*K*E
barB=b*(v-v**-1)/v**2*K.inv()*F
I=sp.eye(n+1)
ss=v**2
check('matrix commutator EF', E*F-F*E-(K-K.inv())/(v-v**-1))
for label,Q,power in [('A',A,-1),('B',B,1),('barA',barA,1),('barB',barB,-1)]:
    check('matrix K weight '+label,K*Q-ss**power*Q*K)
check('matrix AA mixed',A*barA-ss*barA*A-a*(I-K**2))
check('matrix BB mixed',B*barB-ss*barB*B-b*(I-K**-2))
check('matrix A barB mixed',A*barB-ss**-1*barB*A)
check('matrix B barA mixed',B*barA-ss**-1*barA*B)
CC=B*A-ss**-1*A*B
check('matrix triple barA',CC*barA-barA*CC-a*(1-ss**-2)*B)
check('matrix triple barB',CC*barB-barB*CC-b*(ss**-1-ss)*A*K**-2)
for name,X,Y in [('B',B,A),('A',A,B),('barB',barB,barA),('barA',barA,barB)]:
    check('matrix quartic '+name,X**3*Y-(1+ss+ss**-1)*X**2*Y*X
          +(1+ss+ss**-1)*X*Y*X**2-Y*X**3)
check('matrix witness nonzero C',any(sp.cancel(x)!=0 for x in CC))
check('matrix witness nonzero E cubed',any(sp.cancel(x)!=0 for x in E**3))
print('SUMMARY',passed,'PASS, 0 FAIL')
