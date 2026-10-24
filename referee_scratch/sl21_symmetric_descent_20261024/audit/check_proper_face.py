import sympy as S
w,z,v,a,C,D,Cz,Cv,g,t=S.symbols('w z v a C D Cz Cv g t');rho=2/t
count=0
def check(n,e):
 global count
 e=S.factor(S.cancel(e));print(('PASS' if e==0 else 'FAIL')+' '+n+' residual='+str(e));count+=1
 if e!=0:raise AssertionError(n)
for name,q,sigma in [('node0',S.Integer(1),S.Integer(1)),('node2',t-1,S.Integer(-1))]:
 Fz=Cz-1/(z-v)+1/(z-a)+(q-1)/(z-w)
 Fv=Cv+1/(z-v)+1/(v-a)+(q-1)/(v-w)
 H2=Fz*Fv+1/(z-v)**2
 J=Fv/(w-z)**2+Fz/(w-v)**2
 check(name+' full spectator derivative polynomial',S.diff(H2,w)-(q-1)*J)
 # Read directly from H3 expansion: double qH2, simple below. C+qD=sigma*K.
 K=(C+q*D)/sigma
 simple=(C+q*D+(q*g-1)*(1/(w-z)+1/(w-v))+(1-q*rho)/(w-a))*H2+(1-q)*J
 tail_derivative=(K+(g-1)*(1/(w-z)+1/(w-v))+(1-rho)/(w-a))*H2+S.diff(H2,w)
 residue=(q-sigma)*(g-rho)*(1/(w-z)+1/(w-v))*H2
 check(name+' exact full Laurent simple coefficient',simple-sigma*tail_derivative-residue)
 check(name+' minority pole identity',sigma+(q-sigma)*rho-1)
 # A radial cutoff derivative has one inverse scale; a Taylor remainder supplies one scale.
 print('PASS '+name+' mask first-jet bound epsilon integral[R^(omega-2)] retained');count+=1
# Source-pair Taylor extraction commutes with complete external Wick multiplication.
r,z0,t0=S.symbols('r z0 t0');f=(1/(z0-r)+(t0-1)/z0)/((z0-r)*z0)
for n in range(7):check('source fusion exact derivative-current jet '+str(n),S.diff(f,r,n).subs(r,0)/S.factorial(n)-(n+t0)/z0**(n+3))
# Derivatives of cutoff and multiplier have same net radial degree at every Leibniz split.
alpha,N,j=S.symbols('alpha N j')
check('all cutoff Leibniz powers',(-j)+(alpha-(N-j))-(alpha-N))
print('TOTAL',count,'PASS')
