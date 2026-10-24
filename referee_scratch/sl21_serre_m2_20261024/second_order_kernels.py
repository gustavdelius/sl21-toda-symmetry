"""Generate full order-two Wick kernels from charge vectors and cocycles.

The independent unlike-screening exponent is h. No earlier code is imported.
Vertices are (q,m0,m2,has_C), with left momentum q*d+m0*a0+m2*a2
and right momentum m0*a0+m2*a2. Derivative contractions keep t=2/rho.
"""
from itertools import permutations,combinations
import numpy as np
import sympy as S

passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    passed+=int(bool(ok));failed+=int(not bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)

V={'A':(-1,0,0,False),'B':(1,0,0,True),'K':(1,0,1,False),
   'L':(1,1,1,True),'H':(-1,-1,-1,False),
   'Y0':(0,1,0,False),'Y2':(0,0,1,False),'P':(0,1,1,False),'V1':(0,-1,-1,False)}
sources=[.17+.61j,.36-.47j]

def kernel(labels,positions,g,rho,h):
    t=2/rho
    vertices=[V[k] for k in labels]
    C=lambda z:.11+.038*z+.009*z*z
    phi=lambda z:.07*z+.012*z*z
    phi0=lambda z:.11*z+.019*z*z+.003*z**3
    bar0=lambda z:.041*z-.006*z*z
    bar2=lambda z:.037*z+.008*z*z
    normal=0j; power=1.+0j; cocycle=0
    for i,(qi,mi0,mi2,_) in enumerate(vertices):
        zi=positions[i]
        normal+=qi*phi(zi)+mi0*phi0(zi)+mi2*(-t*phi(zi)-phi0(zi))
        normal+=mi0*bar0(np.conj(zi))+mi2*bar2(np.conj(zi))
        for j in range(i+1,len(vertices)):
            qj,mj0,mj2,_=vertices[j]
            full=mi0*mj0+mi2*mj2+(h-1)*(mi0*mj2+mi2*mj0)
            base=rho*qi*qj
            if qi and qj and qi==qj:base+=g-rho
            left=base-qi*(mj0+mj2)-qj*(mi0+mi2)+full
            diff=complex(zi-positions[j])
            power*=diff**left*np.conj(diff)**full
            cocycle+=qi*(mj0+mj2)
    ids=[i for i,v in enumerate(vertices) if v[3]]
    f=[]
    for i in ids:
        val=C(positions[i])
        for j,(qj,mj0,mj2,_) in enumerate(vertices):
            if j!=i:val+=(-qj+mj0+(t-1)*mj2)/(positions[i]-positions[j])
        f.append(val)
    if not ids:W=1
    elif len(ids)==1:W=f[0]
    elif len(ids)==2:W=f[0]*f[1]+1/(positions[ids[0]]-positions[ids[1]])**2
    else:
        W=np.prod(f)+sum(f[k]/(positions[ids[i]]-positions[ids[j]])**2 for i,j,k in ((0,1,2),(0,2,1),(1,2,0)))
    return (-1)**cocycle*np.exp(normal)*power*W

def channel_amplitude(word,coords,g,rho,h,channel):
    t=2/rho;gamma=t/(1-t);base=list(word);B=[i for i,c in enumerate(word) if c=='B'];A=[i for i,c in enumerate(word) if c=='A']
    diagrams=[]
    if channel in ('00','02','22'):
        source_labels={'00':['Y0','Y0'],'02':['Y0','Y2'],'22':['Y2','Y2']}[channel]
        coef=(1 if channel=='02' else .5)/np.pi**2
        diagrams=[(base+source_labels,list(coords)+sources,coef)]
    elif channel in ('K0','K2'):
        for i in B:
            lab=base.copy();lab[i]='K'
            diagrams.append((lab+['Y'+channel[-1]],list(coords)+sources[:1],t/np.pi))
    elif channel=='KK':
        for i,j in combinations(B,2):
            lab=base.copy();lab[i]=lab[j]='K'
            diagrams.append((lab,list(coords),t*t))
    elif channel=='P':diagrams=[(base+['P'],list(coords)+sources[:1],.2**(2*h)/(np.pi*h))]
    elif channel=='L':
        for i in B:
            lab=base.copy();lab[i]='L'
            diagrams.append((lab,list(coords),t*.2**(2*h)/h))
    elif channel=='V1':diagrams=[(base+['V1'],list(coords)+sources[:1],-1/np.pi)]
    elif channel=='H':
        for i in A:
            lab=base.copy();lab[i]='H'
            diagrams.append((lab,list(coords),-gamma))
    return sum(coef*kernel(lab,pos,g,rho,h) for lab,pos,coef in diagrams)

# Check all operator orders against actual ascending-space kernels. Source
# factors and antiholomorphic fractional powers are evaluated, not suppressed.
configs=[(np.sqrt(53.),-np.sqrt(2.),np.sqrt(29.)),
         (2/np.sqrt(5.)+.13,2/np.sqrt(5.),np.sqrt(5.)+.07j),
         (2/np.sqrt(5.),2/np.sqrt(5.),np.sqrt(5.))]
xs=np.array([-.83,-.24,.31,.89]);perms=list(permutations(range(4)))
for setting,(g,rho,h) in enumerate(configs):
    u=np.exp(-1j*np.pi*g);s=np.exp(-1j*np.pi*rho)
    co=np.array([1,-1-s-1/s,1+s+1/s,-1])
    defect=(s-1)*(s+1)*(s-u)*(s*s-u)*(u+1)/s**4
    for major in ('B','A'):
        minor='A' if major=='B' else 'B'
        words=[major*(3-k)+minor+major*k for k in range(4)]
        channels=['00','02','22','K0','K2','KK','P','L','V1','H']
        if major=='A':channels.remove('KK')
        for channel in channels:
            canonical=np.array([channel_amplitude(word,xs,g,rho,h,channel) for word in words])
            scale=max(abs(canonical))
            sums=np.zeros(4,complex);worst=0
            for col,word in enumerate(words):
                for order in perms:
                    rank={i:k for k,i in enumerate(order)}
                    coords=[xs[rank[i]] for i in range(4)]
                    row=words.index(''.join(word[i] for i in order))
                    factor=1.+0j
                    for i in range(4):
                        for j in range(i+1,4):
                            if rank[i]>rank[j]:factor*=u if word[i]==word[j] else 1/s
                    actual=channel_amplitude(word,coords,g,rho,h,channel)
                    target=factor*canonical[row]
                    worst=max(worst,abs(actual-target)/(scale+abs(target)+1e-200))
                    sums[col]+=actual
            check(f'full order-two Wick exchange {major}3 {channel}, setting {setting}',worst<5e-12,
                  f'error={worst:.3g}')
            target=defect*(canonical[2]-canonical[1])
            error=abs(np.dot(co,sums)-target)/(scale+abs(target)+1e-200)
            check(f'quartic chamber formula {major}3 {channel}, setting {setting}',error<2e-10,
                  f'error={error:.3g}')

# Reconstruct the crucial first-order scalar derivative from Wick's theorem.
a,w,t=S.symbols('a w t',nonzero=True);rho=2/t;beta=t-1;eta=t-2
z=[S.Rational(-7,8),S.Rational(1,6),S.Rational(11,9)]
f=[-sum(1/(zi-zj) for j,zj in enumerate(z) if i!=j)+1/(zi-a)+beta/(zi-w) for i,zi in enumerate(z)]
W=S.prod(f)+sum(f[k]/(z[i]-z[j])**2 for i,j,k in ((0,1,2),(0,2,1),(1,2,0)))
G=(a-w)**2/S.prod((w-zi)**2 for zi in z)
rhs=beta*eta/rho*(S.diff(G,a)+rho*sum(1/(zi-a) for zi in z)*G)
num=S.together(W*(a-w)/S.prod(zi-w for zi in z)-rhs).as_numer_denom()[0]
check('fresh first-order marginal derivative audit',S.expand(num)==0)
print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
