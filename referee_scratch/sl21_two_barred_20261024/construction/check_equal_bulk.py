import sympy as s
k,r,q,g,t=s.symbols('k r q g t')
A=1/r**2;B=s.I*k/r-1/r**2
P=2*q*(q-1)*A*B+q*q*(A*A+B*B)
assert s.cancel(P-(q*q*(s.I*k)**2/r**2-2*q*s.I*k/r**3+2*q/r**4))==0
print('PASS complete equal-source current-integral residue polynomial')
coef=q*q*g*(g-1)-2*q*g+2*q
assert s.factor(coef)==q*(g-1)*(g*q-2)
print('PASS equal-source current-integral coefficient q(g-1)(qg-2)')
assert s.cancel((coef-(1-g)*(2-g)).subs(q,2/g-1))==0
print('PASS physical equal-source coefficient (1-g)(2-g)')
print('PASS total 3')
