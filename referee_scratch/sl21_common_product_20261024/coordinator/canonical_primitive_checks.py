"""Independent scalar tail primitives and canonical/action pair comparison."""
import sympy as s
t,a,w=s.symbols('t a w', nonzero=True)
rho=2/t; eta=t-2;gamma=t/(1-t)
z=[s.Integer(-2),s.Integer(1),s.Integer(4)]

def ck(name,v):
 assert s.factor(v)==0,(name,s.factor(v))
 print('PASS '+name)

for i in range(3):
 j,k=[j for j in range(3) if j!=i]
 Li=1/(z[j]-a)+eta/(z[j]-z[i])-1/(z[j]-z[k])
 Lj=1/(z[k]-a)+eta/(z[k]-z[i])-1/(z[k]-z[j])
 poly=Li*Lj+1/(z[j]-z[k])**2
 extra=(z[i]-a)/((z[i]-z[j])*(z[i]-z[k]))
 Si=(a-z[i])**2/((z[i]-z[j])**2*(z[i]-z[k])**2)
 logRa=-rho*sum(1/(a-zi) for zi in z)
 ck('individual K2 tail minority derivative '+str(i),extra*poly+eta/rho*(s.diff(Si,a)+logRa*Si))

for i in range(3):
 Si=(w-z[i])**2/s.prod((z[i]-z[j])**2 for j in range(3) if j!=i)
 derivative=sum(1/(w-zi) for zi in z)-t/(w-z[i])
 logRw=-rho*sum(1/(w-zi) for zi in z)
 ck('individual F_A tail minority derivative '+str(i),Si*derivative+(s.diff(Si,w)+logRw*Si)/rho)

# Direct two-tail scalar contractions: each B/K2 scalar factor is
# (t-2)/(z-K), the additional source exponential costs 1/(K-z).
# With its -t coupling coefficient each tail contributes +t(t-2)/r².
canonical_BB=2*t*(t-2)
action_BB=-2*(t-1)*(t-2)
ck('BB canonical plus action pair coefficient',canonical_BB+action_BB-2*(t-2))
canonical_AA=-2*gamma
action_AA=-2
ck('AA canonical plus action pair coefficient',canonical_AA+action_AA+2/(1-t))
print('INFO canonical/action BB ratio',s.factor(canonical_BB/(canonical_BB+action_BB)))
print('INFO canonical/action AA ratio',s.factor(canonical_AA/(canonical_AA+action_AA)))
print('INFO The two products cannot be identified termwise; an interacting map is additional data.')
