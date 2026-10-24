import sympy as s
h,t,q,qb,e=s.symbols('h t q qb epsilon',nonzero=True)
eta=h-t
R=q**eta*qb**eta
for i,j in [(1,0),(2,0),(0,1),(0,2),(1,1)]:
 # Derive falling-factorial coefficient without branch-sensitive simplification.
 actual=s.prod(eta-k for k in range(i))*s.prod(eta-k for k in range(j))
 quotient=s.cancel(actual/eta)
 assert s.denom(quotient)==1
 print('PASS differentiated regulator divisible by h-t:',i,j, 'quotient=',s.expand(quotient))
checks={
 'action CT transport':(-1/h+1/t)-eta/(h*t),
 'tail CT transport':t/h-1+eta/h,
 'endpoint scalar CT transport':(-1/(h-1)+1/(t-1))-eta/((h-1)*(t-1)),
 'endpoint descendant CT transport':-1/h+1/t-eta/(h*t),
}
for name,x in checks.items():
 assert s.factor(x)==0,name
 print('PASS '+name)
# Formal second-derivative Leibniz expansion with ANY smooth cutoff (including approximated masks).
Q=s.Function('Q')(q);M=s.Function('M')(q)
x=s.diff(R*M*Q,q,2)-R*s.diff(M*Q,q,2)
expected=2*s.diff(R,q)*s.diff(M*Q,q)+s.diff(R,q,2)*M*Q
assert s.simplify(x-expected)==0
print('PASS arbitrary-mask second derivative defect is only differentiated regulator')
print('Scope: distribution product rule in initial integrable domain; no assertion of full associativity.')
