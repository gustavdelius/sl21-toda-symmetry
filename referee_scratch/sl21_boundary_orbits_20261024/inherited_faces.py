"""Fresh symbolic Wick audit of the inherited full marginal channel."""
import sympy as S

z=list(S.symbols('z1:4'))
a,w,t=S.symbols('a w t',nonzero=True)
rho=2/t
n=0
def check(name,expr):
    global n
    result=S.factor(expr)
    if result!=0:
        print('FAIL',name,result)
        raise AssertionError(name)
    n+=1
    print('PASS',name)

def scalar(fs,positions):
    if len(fs)==3:
        return S.prod(fs)+sum(fs[i]/(positions[(i+1)%3]-positions[(i+2)%3])**2 for i in range(3))
    return S.prod(fs)+1/(positions[0]-positions[1])**2

f0=[-sum(1/(x-y) for y in z if y!=x)+1/(x-a) for x in z]
check('BBBA source-free scalar cancels',scalar(f0,z))
linear=sum(f0[i]*f0[j]+1/(z[i]-z[j])**2 for i in range(3) for j in range(i+1,3))
check('BBBA source-free first common C jet cancels',linear)
f2=[f+(t-1)/(x-w) for f,x in zip(f0,z)]
f1=[f+1/(x-w) for f,x in zip(f0,z)]
check('BBBA node-zero scalar cancels',scalar(f1,z))
G=(a-w)**2/S.prod((w-x)**2 for x in z)
source=(a-w)/S.prod(x-w for x in z)*scalar(f2,z)
logR=rho*sum(1/(x-a) for x in z)
check('BBBA node-two scalar is minority-coordinate derivative',
      source-(t-1)*(t-2)/rho*(S.diff(G,a)+logR*G))

for i in range(3):
    ai=z[i]
    Si=(a-ai)**2/S.prod((ai-x)**2 for x in z if x!=ai)
    others=[j for j in range(3) if j!=i]
    # Direct scalar contractions for the two remaining B derivatives with
    # a K tail at ai, the other B, and the A current at a.
    fk=[]
    for j in others:
        x=z[j]
        other=z[others[1] if j==others[0] else others[0]]
        fk.append(-1/(x-other)+(t-2)/(x-ai)+1/(x-a))
    vertex=(a-ai)/S.prod(x-ai for x in z if x!=ai)
    # Current-first cocycles contribute -1 for each remaining current to K;
    # the three signs combine to -1, cancelling the tail's -t coefficient.
    actual=t*vertex*scalar(fk,[z[j] for j in others])
    check(f'BBBA K tail at slot {i} minority derivative',
          actual-t*(t-2)/rho*(S.diff(Si,a)+logR*Si))

zz=z[:2]
ff=[-1/(x-y)+1/(x-a) for x,y in [zz,zz[::-1]]]
check('proper BBA scalar cancels',scalar(ff,zz))
print(f'{n} PASS, 0 FAIL')
