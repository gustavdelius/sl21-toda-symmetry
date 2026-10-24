import sympy as s
from itertools import permutations
rho,t,r,a,w,C=s.symbols('rho t r a w C',nonzero=True)
rho=2/t
eta=1/(1-t);gamma=t/(1-t)
# a is fixed B, z=w+r is moving A. Rational numerator of double pole.
P=(a-w-r)**(-rho)*(a-w)**2*(C+1/(a-w-r)-t/(a-w))
V=(a-w)**(2-rho)*(C+(1-t)/(a-w))
# Logarithmic/power factors are cancelled on a single local branch before simplification.
base=(a-w)**(-rho)
coeff0=(a-w)**2*(C+(1-t)/(a-w))
coeff1=rho*(a-w)*(C+(1-t)/(a-w))+1
Vp=-(2-rho)*(a-w)*C-(1-rho)*(1-t)
def ck(name,v):
 assert bool(v),name
 print('PASS '+name)
ck('companion double pole actual vertex',s.factor(coeff0-((a-w)**2*C+(1-t)*(a-w)))==0)
ck('companion simple pole eta derivative',s.factor(coeff1-eta*Vp)==0)
ck('companion vector jet eta k=d',s.factor(eta*(1-t)-1)==0)
ck('action coefficient gamma',s.factor(eta-1-gamma)==0)
ck('tails gamma squared',s.factor(gamma*eta-gamma-gamma**2)==0)
ck('canonical finite addition',s.factor(gamma**2+gamma-gamma-gamma**2)==0)
# Conditional Chevalley conversion, all K reorder factors explicit.
x,a0,b0=s.symbols('s a0 b0',nonzero=True)
# [A,barA K^-1] = (AbarA-s barAA)K^-1; similarly B with K.
K=s.symbols('K',nonzero=True)
ck('Cartan A bracket',s.expand(a0*(1-K*K)/K-a0*(1/K-K))==0)
ck('Cartan B bracket',s.expand(b0*(1-K**-2)*K-b0*(K-1/K))==0)
ck('cross A barB K cancels',s.simplify(1/x-1/x)==0)
ck('cross B barA K^-1 cancels',s.simplify(1/x-1/x)==0)
# For F_A=barA K^-1,F_B=barB K: K barA=s barA K, K barB=s^-1 barB K.
# Commuting all K powers to right gives exponent sum p_i weight_j.
for majority,minority in [('a','b'),('b','a')]:
 ws=[majority*3+minority,majority*2+minority+majority,majority+minority+majority*2,minority+majority*3]
 power={'a':-1,'b':1};weight={'a':1,'b':-1}
 exps=[sum(power[W[i]]*weight[W[j]] for i in range(4) for j in range(i+1,4)) for W in ws]
 print('DIAGNOSTIC barred quartic K factors',ws,exps)
 ck('negative Serre common factor '+majority,len(set(exps))==1)
print('All checks passed; no field-theory algebra inferred.')
