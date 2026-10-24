"""Direct three-dimensional endpoint integral after exact source/A integration.

Retains the contour-height terms and the finite circular holes at inner
currents. No earlier script or formula implementation is imported.
"""
import numpy as np
from scipy.stats import qmc
from pathlib import Path

rho=np.sqrt(3.)
t=2/rho
beta=t-1
s=np.exp(-1j*np.pi*rho)
cs=1+s+1/s
coef=np.array([1,-cs,cs,-1])
height_slots=np.array([1.5,.5,-.5,-1.5])
c=np.sqrt(5.)

def primitive(x,ep,k,end,holes=True,heights=True):
    apos=3-k
    slots=[j for j in range(4) if j!=apos]
    z=x+1j*c*ep*height_slots[slots]
    a=end+1j*c*ep*height_slots[apos]
    R=np.ones(len(x),complex)
    for i in range(3):
        R*=np.power(z[:,i]-a,-rho)
        for j in range(i+1,3):
            R*=np.power(z[:,i]-z[:,j],rho)
    Sh=np.zeros(len(x),complex)
    H=np.zeros(len(x),complex)
    T=np.zeros(len(x),complex)
    for i in range(3):
        others=[j for j in range(3) if j!=i]
        dif=[z[:,i]-z[:,j] for j in others]
        D=1/(dif[0]**2*dif[1]**2)
        l=1/dif[0]+1/dif[1]
        m=1/dif[0]**2+1/dif[1]**2
        v=a-z[:,i]
        Sh+=v*v*D
        Ri=-2*D*(v+v*v*l)
        H+=z[:,i].imag*Ri
        if slots[i]!=0:
            T+=D*(1+4*v*l+v*v*(2*l*l+m))
    return s**k*(t-2)/rho*R*(Sh-(2j*beta*H if heights else 0)+(beta*ep*ep*T if holes else 0))

def run(seed,power=16):
    x=2*qmc.Sobol(d=3,scramble=True,seed=seed).random_base2(power)-1
    ans=[]
    for ep in (.04,.02,.01,.005,.0025):
        values=[]
        for k in range(4):
            values.append(8*np.mean(primitive(x,ep,k,2)-primitive(x,ep,k,-2)))
        values=np.array(values)
        total=np.dot(coef,values)
        ans.append((ep,total))
        print(f'DATA seed={seed} epsilon={ep:.6g} sum={total.real:.12g}{total.imag:+.12g}j abs={abs(total):.8g}',flush=True)
    return ans

data=[run(1729),run(9323)]
passed=failed=0
def check(name,ok,detail=''):
    global passed,failed
    if ok:passed+=1
    else:failed+=1
    print(('PASS ' if ok else 'FAIL ')+name+(': '+detail if detail else ''),flush=True)

check('endpoint phase polynomial at s',abs(sum(coef[k]*s**k for k in range(4)))<2e-14)
check('endpoint phase polynomial at inverse s',abs(sum(coef[k]*s**(-k) for k in range(4)))<2e-14)
for seed,rows in zip((1729,9323),data):
    check(f'finite-height/hole boundary sum decreases, seed {seed}',abs(rows[-1][1])<abs(rows[0][1])/3)
    check(f'last two cutoff values decrease, seed {seed}',abs(rows[-1][1])<abs(rows[-2][1]))
spread=abs(data[0][-1][1]-data[1][-1][1])
check('independent QMC agreement at smallest cutoff',spread<.15*max(abs(data[0][-1][1]),abs(data[1][-1][1])),f'spread={spread:.6g}')
print(f'RESULT {passed} PASS, {failed} FAIL',flush=True)
raise SystemExit(bool(failed))
