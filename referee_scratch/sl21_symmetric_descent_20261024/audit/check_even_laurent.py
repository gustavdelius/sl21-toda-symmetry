import sympy as S
t,g,D=S.symbols('t g D');rho=2/t;sig=1/(1-t)
checks=[('even Laurent normal field',-D-sig*(t-1)*D),('even Laurent majority coefficient',g-sig*(g-2)-(1-sig)*(g-rho)),('even Laurent minority coefficient',-rho-sig*(2-rho)),('even Laurent derivative spectator',1-sig*(1-t))]
for n,e in checks:
 e=S.factor(e);print(('PASS' if e==0 else 'FAIL')+' '+n+' residual='+str(e));assert e==0
print('TOTAL 4 PASS')
